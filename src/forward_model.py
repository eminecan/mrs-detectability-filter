"""
forward_model.py
================
Hyperpolarized 13C pyruvate -> {lactate, alanine, bicarbonate} exchange model,
plus the detectability threshold that turns a transcript fold-change into a
"can I actually see this in vivo?" verdict.

This is the bridge from Perturb-seq (transcript log-fold-change) to MRS (flux).
The physics here is generic and citable; the biology (which gene sets a rate
constant, and by how much) lives in config/gene_sets.yaml and is authored by
the researcher.

Standard unidirectional multi-product model for HP substrates (reverse fluxes
neglected because product pools are small over the acquisition window):

    dP/dt  = -(kPL + kPA + kPB + R1p) * P
    dL/dt  =  kPL * P - R1l * L
    dA/dt  =  kPA * P - R1a * A
    dB/dt  =  kPB * P - R1b * B

with R1x = 1 / T1x. Default T1 values are order-of-magnitude in vivo figures at
3T; REPLACE with the values from your own acquisitions before quoting numbers.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp

# NumPy 2.x renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# ---------------------------------------------------------------------------
# Default relaxation / kinetic parameters.  TODO(you): set from your own data.
# ---------------------------------------------------------------------------
@dataclass
class HPParams:
    # Longitudinal relaxation times (s) — in vivo, 3T, approximate.
    T1_pyr: float = 30.0
    T1_lac: float = 25.0
    T1_ala: float = 25.0
    T1_bic: float = 20.0
    # Baseline apparent rate constants (1/s) for the *unperturbed* cell state.
    # These are the numbers your model perturbs. TODO(you): anchor to literature
    # or to your Sci Rep 2020 activated-T-cell measurements.
    kPL: float = 0.030
    kPA: float = 0.005
    kPB: float = 0.004
    # Acquisition
    t_end: float = 60.0      # s
    dt: float = 0.5          # s
    P0: float = 1.0          # normalized initial polarized pyruvate

    @property
    def R1(self):
        return dict(pyr=1/self.T1_pyr, lac=1/self.T1_lac,
                    ala=1/self.T1_ala, bic=1/self.T1_bic)


def simulate(params: HPParams, kPL=None, kPA=None, kPB=None):
    """Integrate the exchange ODEs. Returns (t, dict of concentration curves)."""
    kPL = params.kPL if kPL is None else kPL
    kPA = params.kPA if kPA is None else kPA
    kPB = params.kPB if kPB is None else kPB
    R1 = params.R1

    def rhs(t, y):
        P, L, A, B = y
        dP = -(kPL + kPA + kPB + R1['pyr']) * P
        dL = kPL * P - R1['lac'] * L
        dA = kPA * P - R1['ala'] * A
        dB = kPB * P - R1['bic'] * B
        return [dP, dL, dA, dB]

    t = np.arange(0, params.t_end + params.dt, params.dt)
    sol = solve_ivp(rhs, (0, params.t_end), [params.P0, 0, 0, 0],
                    t_eval=t, method='LSODA', rtol=1e-8, atol=1e-10)
    P, L, A, B = sol.y
    return t, dict(pyr=P, lac=L, ala=A, bic=B)


def lac_pyr_ratio_auc(params: HPParams, kPL=None):
    """A common in vivo observable: AUC(lactate)/AUC(pyruvate)."""
    t, c = simulate(params, kPL=kPL)
    return _trapz(c['lac'], t) / _trapz(c['pyr'], t)


# ---------------------------------------------------------------------------
# Transcript fold-change -> predicted change in a rate constant.
# ---------------------------------------------------------------------------
def scale_rate_by_node(baseline_k: float, log2fc: float,
                       sign: int = 1, coupling: float = 1.0) -> float:
    """
    Map a flux-determining node's log2 fold-change to a new rate constant.

    sign (+1/-1) encodes the direction of the node's control over the rate
    constant. effective = sign * log2fc, so a knockdown (log2fc < 0) of a -1
    node (inhibitor / reverse-favoring enzyme, e.g. PDK, GLUL, LDHB) RAISES the
    rate constant, as it should.

    coupling in [0,1] is a deliberately conservative factor: transcript change
    rarely translates 1:1 to enzyme activity / flux. Report the number AND the
    coupling you assumed. Do not pretend coupling = 1 is measured truth. Use a
    lower coupling for post-translationally regulated arms (e.g. PDK on kPB).
    """
    effective = sign * log2fc
    return baseline_k * (2.0 ** effective) ** coupling


# ---------------------------------------------------------------------------
# THE detectability verdict.  This is the line the project is built around.
# ---------------------------------------------------------------------------
@dataclass
class DetectabilityResult:
    rate_constant: str
    k_baseline: float
    k_perturbed: float
    frac_change: float
    detectable: bool
    threshold: float
    note: str = ""


def detectability(baseline_k: float, perturbed_k: float, rate_name: str,
                  threshold: float = 0.20) -> DetectabilityResult:
    """
    In vivo HP-13C differences below ~15-25% in a rate constant are not
    reliably resolvable. Default 0.20. A candidate below threshold is
    'mechanistically real, not PD-trackable' — and you SAY SO on camera.
    """
    frac = abs(perturbed_k - baseline_k) / baseline_k
    ok = frac >= threshold
    note = ("trackable by HP-13C MRS" if ok
            else f"predicted delta {frac:.0%} < {threshold:.0%} floor — "
                 "real biology, not PD-trackable")
    return DetectabilityResult(rate_name, baseline_k, perturbed_k,
                               frac, ok, threshold, note)


if __name__ == "__main__":
    p = HPParams()

    # Example 1: LDHA knockdown (sign +1, forward driver), log2fc=-1, coupling 0.7
    kPL_new = scale_rate_by_node(p.kPL, log2fc=-1.0, sign=+1, coupling=0.7)
    r1 = detectability(p.kPL, kPL_new, "kPL")
    print("LDHA knockdown (sign +1):")
    print(f"  kPL {r1.k_baseline:.4f} -> {r1.k_perturbed:.4f}  "
          f"AUC(Lac/Pyr) {lac_pyr_ratio_auc(p):.3f} -> "
          f"{lac_pyr_ratio_auc(p, kPL=kPL_new):.3f}")
    print(f"  {r1.note}\n")

    # Example 2: PDK knockdown (sign -1, inhibitor) should RAISE kPB.
    kPB_new = scale_rate_by_node(p.kPB, log2fc=-1.0, sign=-1, coupling=0.4)
    r2 = detectability(p.kPB, kPB_new, "kPB")
    print("PDK knockdown (sign -1, low coupling for post-translational arm):")
    print(f"  kPB {r2.k_baseline:.4f} -> {r2.k_perturbed:.4f}   ({r2.note})")
