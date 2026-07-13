"""Called by generate_mosqito_zwtv_crosscheck.jl — not directly.

Reads {name, fs, field_type, signal} cases from JSON (argv[1]), runs
MoSQITo's loudness_zwtv @ d990c33f94f1 (Apache-2.0) both as the public
end-to-end API and, per case, by calling its private per-stage functions
directly (band levels -> core loudness -> nonlinear decay -> spreading
-> temporal weighting -> decimation), and writes a Julia fixtures file
(argv[2]) with the per-case summary plus full stage intermediates for
one designated case (argv[3], a case name).

Optional argv[4]/argv[5]/argv[6]: absolute path to a *local, uncommitted*
ISO 532-1 Annex B time-varying .wav file (from a MoSQITo checkout), a
fixture name, and a provenance string. If given, this end-to-end result
is ALSO written to argv[2] but with ONLY the final derived numbers (N_t,
time_axis, N5, N10) -- deliberately excluding band_levels and
specific_loudness, which (unlike a handful of final loudness numbers)
retain enough structure to reconstitute the copyrighted recording. See
.superpowers/sdd/zwtv-pins.md for the licensing reasoning. The .wav
itself is never read by, nor shipped with, this repository's committed
sources -- it must be supplied from a local MoSQITo checkout at
generation time.

Stage boundary transcribed from loudness_zwtv.py:111-131 (module docstring
Notes therein cite ISO 532-1:2017 clause 6 "Method for time-varying
sounds"; _third_octave_levels.py:17 cites "ISO 532-1 section 6.3"):
    band_levels        = _third_octave_levels(sig, fs)[0]          # 28 x nblocks @ 2000 Hz
    nm_pre_decay        = _main_loudness(band_levels, field_type)   # 21 x nblocks
    nm_post_decay       = _nl_loudness(nm_pre_decay)                # 21 x nblocks (KERNEL: nonlinear decay)
    N_pre_weight, Nsp   = _calc_slopes(nm_post_decay)                # nblocks, / 240 x nblocks (KERNEL: spreading)
    N_pre_decimation    = _temporal_weighting(N_pre_weight)         # nblocks (KERNEL: temporal weighting)
    N_final             = N_pre_decimation[::4]                     # KERNEL: dec_factor=4 (loudness_zwtv.py:128-131)
    Nsp_final           = Nsp[:, ::4]                                # NOTE: specific loudness is decimated but
                                                                      # NEVER temporally weighted -- only total N is.

For every case we also call the public `loudness_zwtv` and assert its N
matches our per-stage reconstruction (N_final) to machine precision --
a mismatch would mean this driver mis-transcribed the private call
sequence, not a MoSQITo bug.
"""
import hashlib
import json
import sys

import numpy as np
from mosqito.sq_metrics import loudness_zwtv
from mosqito.sq_metrics.loudness.loudness_zwst._calc_slopes import _calc_slopes
from mosqito.sq_metrics.loudness.loudness_zwst._main_loudness import _main_loudness
from mosqito.sq_metrics.loudness.loudness_zwtv._nonlinear_decay import _nl_loudness
from mosqito.sq_metrics.loudness.loudness_zwtv._temporal_weighting import (
    _temporal_weighting,
)
from mosqito.sq_metrics.loudness.loudness_zwtv._third_octave_levels import (
    _third_octave_levels,
)


def jl_vec(v):
    return "[" + ", ".join(repr(float(x)) for x in np.asarray(v).ravel()) + "]"


def jl_mat(m):
    # Julia matrix literal, row-major source (numpy) -> Julia's row-major
    # literal syntax `[a b c; d e f]` (rows separated by ';', cols by ' ').
    rows = ["  " + " ".join(repr(float(x)) for x in row) for row in np.asarray(m)]
    return "[\n" + ";\n".join(rows) + "\n]"


def sha256_of(arr):
    return hashlib.sha256(np.asarray(arr, dtype=np.float64).tobytes()).hexdigest()


cases = json.load(open(sys.argv[1]))
out_path = sys.argv[2]
full_dump_case = sys.argv[3]

case_entries = []
stage_block = None

for c in cases:
    name = c["name"]
    fs = float(c["fs"])
    field_type = c["field_type"]
    sig = np.array(c["signal"], dtype=np.float64)

    # --- Per-stage reconstruction (private functions) ---
    band_levels, band_time_axis, _freq = _third_octave_levels(sig, fs)
    nm_pre_decay = _main_loudness(band_levels, field_type)
    nm_post_decay = _nl_loudness(nm_pre_decay)
    N_pre_weight, N_specific_pre_weight = _calc_slopes(nm_post_decay)
    N_pre_decimation = _temporal_weighting(N_pre_weight)

    dec_factor = 4
    N_final = N_pre_decimation[::dec_factor]
    N_specific_final = N_specific_pre_weight[:, ::dec_factor]
    time_axis_final = band_time_axis[::dec_factor]

    # --- Public API cross-check: must match our reconstruction exactly ---
    N_ref, N_spec_ref, _bark_ref, time_ref = loudness_zwtv(sig, fs, field_type)
    assert np.allclose(N_ref, N_final, rtol=0, atol=0) or np.array_equal(
        N_ref, N_final
    ), f"{name}: public API N diverges from private-stage reconstruction"
    assert np.array_equal(N_spec_ref, N_specific_final), (
        f"{name}: public API N_specific diverges from private-stage reconstruction"
    )
    assert np.array_equal(time_ref, time_axis_final), (
        f"{name}: public API time_axis diverges from private-stage reconstruction"
    )

    n5 = float(np.percentile(N_final, 95, method="linear"))
    n10 = float(np.percentile(N_final, 90, method="linear"))

    case_entries.append(
        f"""    (
        name = "{name}",
        fs = {fs!r},
        field_type = :{field_type},
        band_levels_sha256 = "{sha256_of(band_levels)}",
        band_levels = {jl_mat(band_levels)},
        N_t = {jl_vec(N_final)},
        time_axis = {jl_vec(time_axis_final)},
        N5_iso = {n5!r},
        N10_iso = {n10!r},
    ),"""
    )

    if name == full_dump_case:
        # Only the genuinely NEW intermediates go here -- band_levels/N_t/
        # time_axis for this case already live in ZWTV_KERNEL_FIXTURE_CASES
        # (matched by `name`); duplicating them would just bloat the file.
        # specific_loudness (240 x nblocks) is intentionally NOT vendored:
        # not required by the plan's fixture shape and, at any usable
        # duration, dwarfs everything else here -- Task 3's kernel-level
        # conformance pass (Annex B band-level inputs) is the place for
        # that, generated fresh from the vendored band_levels if needed.
        stage_block = f"""
# Full per-stage intermediates for case "{name}" (see that case in
# ZWTV_KERNEL_FIXTURE_CASES for its band_levels input and final N_t/
# time_axis) -- for unit-testing each kernel stage in isolation (Task 2).
# Stage order and citations as in the module docstring above.
const ZWTV_STAGE_INTERMEDIATES = (
    name = "{name}",
    nm_post_decay = {jl_mat(nm_post_decay)},           # _nonlinear_decay output: 21 x nblocks
    N_pre_decimation = {jl_vec(N_pre_decimation)},     # post-temporal-weighting, pre-decimation N(t)
)
"""

assert stage_block is not None, f"full_dump_case {full_dump_case!r} not found in cases"

# --- Percentile-definition pin: ISO/DIN (Hyndman-Fan type 7, numpy/Julia
# default "linear") vs MATLAB prctile (Hazen, type 5, numpy method="hazen").
# Fixed vector, independent of any synthesized signal, for a reproducible,
# generator-independent pin.
PIN_VECTOR = [
    0.12, 0.5, 1.3, 2.75, 3.0, 3.4, 4.125, 4.9, 5.0, 5.5,
    6.2, 6.75, 7.1, 7.8, 8.25, 8.9, 9.4, 9.75, 10.1, 10.6,
]
pin = np.array(PIN_VECTOR)
n5_type7 = float(np.percentile(pin, 95, method="linear"))
n10_type7 = float(np.percentile(pin, 90, method="linear"))
n5_type5 = float(np.percentile(pin, 95, method="hazen"))
n10_type5 = float(np.percentile(pin, 90, method="hazen"))

pin_block = f"""
# Percentile-definition pin (Task 2 consumes this). ISO 532-1/DIN 45631/A1
# define N5 as the loudness reached or exceeded 5 % of the time (95th
# percentile of N(t)), N10 analogously at 90 %. Julia's `quantile` default
# and numpy's `percentile` default are both Hyndman-Fan type 7 ("linear" in
# numpy's `method` kwarg) -- no divergence there. MATLAB's `prctile` (used
# by SQAT, the eventual PA/N5 crosscheck) instead follows a Hazen-style
# scheme, numpy `method="hazen"` (Hyndman-Fan type 5) approximates it. This
# fixed vector records the measured gap for later attribution; it is NOT
# derived from any copyrighted signal.
const ZWTV_PERCENTILE_PIN_VECTOR = {jl_vec(pin)}
const ZWTV_PERCENTILE_PIN = (
    n5_type7 = {n5_type7!r},   # ISO/DIN definition we ship (95th pct, type 7 / "linear")
    n10_type7 = {n10_type7!r}, # ISO/DIN definition we ship (90th pct, type 7 / "linear")
    n5_type5 = {n5_type5!r},   # MATLAB-prctile-style proxy (Hazen / type 5), for gap attribution only
    n10_type5 = {n10_type5!r},
)
"""

header = '''# AUTOGENERATED by scripts/generate_mosqito_zwtv_crosscheck.jl -- do not edit.
# (case, fs, field_type, band_levels [28 x nblocks @ 2000 Hz], N_t [sone,
# 2 ms axis], time_axis [s], N5_iso/N10_iso [sone, ISO/DIN type-7
# percentile]) computed by MoSQITo loudness_zwtv @ d990c33f94f1
# (Apache-2.0) via its private per-stage functions, on signals synthesized
# by this package's own generator script -- any disagreement is a
# transcription bug. Every case's public-API N/N_specific/time_axis was
# asserted equal to the private-stage reconstruction at generation time
# (see crosscheck_zwtv.py).
#
# NOTE on N5_iso/N10_iso float precision: these were computed by numpy
# percentile(method="linear") on the numpy-side N(t); Julia
# `Statistics.quantile` on the same vendored N_t can differ by 1 ULP
# (measured: exactly 1 ULP on N5_iso for all three cases, 0 ULP on
# N10_iso -- summation-order difference in the interpolation arithmetic,
# not a definition mismatch). Tests must NOT assert bit-exact equality on
# these two fields; use rtol ~ 4*eps() or compare against Julia-side
# quantile of N_t. The dedicated ZWTV_PERCENTILE_PIN below IS bit-exact
# between numpy and Julia (verified) -- use it for the definition pin.
const ZWTV_KERNEL_FIXTURE_CASES = [
'''

annexb_block = ""
if len(sys.argv) > 4:
    wav_path, annexb_name, provenance = sys.argv[4], sys.argv[5], sys.argv[6]
    from mosqito.utils import load as mosqito_load

    # Calibration factor transcribed from MoSQITo's own Annex B validation
    # driver (validations/sq_metrics/loudness_zwtv/validation_loudness_zwtv.py:191)
    # so our numbers reproduce the same reference-comparable results.
    sig, fs = mosqito_load(wav_path, wav_calib=2 * 2**0.5)
    N_annexb, N_spec_annexb, _bark_annexb, time_annexb = loudness_zwtv(
        sig, float(fs), "free"
    )

    # Same public-API-vs-private-stage consistency assertion as the
    # synthetic-case loop above -- the Annex B case is not exempt.
    bl_b, bt_b, _freq_b = _third_octave_levels(sig, float(fs))
    nm_b = _nl_loudness(_main_loudness(bl_b, "free"))
    N_pw_b, Nsp_b = _calc_slopes(nm_b)
    N_rec_b = _temporal_weighting(N_pw_b)[::4]
    assert np.array_equal(N_annexb, N_rec_b), (
        f"{annexb_name}: public API N diverges from private-stage reconstruction"
    )
    assert np.array_equal(N_spec_annexb, Nsp_b[:, ::4]), (
        f"{annexb_name}: public API N_specific diverges from private-stage reconstruction"
    )
    assert np.array_equal(time_annexb, bt_b[::4]), (
        f"{annexb_name}: public API time_axis diverges from private-stage reconstruction"
    )

    n5_annexb = float(np.percentile(N_annexb, 95, method="linear"))
    n10_annexb = float(np.percentile(N_annexb, 90, method="linear"))
    # Coarse decimation (~25x) of the already-2-ms N(t)/time_axis for a
    # lightweight shape/regression check -- N5/N10 above are computed on
    # the FULL series, so the percentile numbers are faithful; only the
    # illustrative curve itself is thinned, keeping this "scratch" case
    # small (full-resolution N(t) for a 10s+ signal would dominate the
    # file for a fixture the plan calls optional/scratch).
    coarse = 25
    annexb_block = f"""
# ISO 532-1 Annex B time-varying conformance case, DERIVED NUMBERS ONLY.
# Provenance: {provenance}
# The underlying .wav is ISO-copyrighted material MoSQITo redistributes
# under its own Apache-2.0 grant, which cannot confer rights ISO has not
# granted -- so the recording itself is never read by, nor shipped with,
# this repository. Only the final loudness curve (coarsened {coarse}x; N5/
# N10 below are exact, computed on the full series) and its ISO/DIN
# percentiles (a highly-compressed perceptual summary, not a
# reconstruction of the recording) are vendored here, exactly as
# test/fixtures/test_signal_1_nspec.csv already vendors MoSQITo-computed
# Annex B *results* for Method 1. band_levels/specific_loudness are
# deliberately withheld for this case (see zwtv-pins.md).
const ZWTV_ANNEXB_DERIVED = (
    name = "{annexb_name}",
    N_t_coarse = {jl_vec(N_annexb[::coarse])},
    time_axis_coarse = {jl_vec(time_annexb[::coarse])},
    coarse_factor = {coarse},
    n_full = {len(N_annexb)},
    N5_iso = {n5_annexb!r},
    N10_iso = {n10_annexb!r},
)
"""

with open(out_path, "w") as f:
    f.write(header)
    f.write("\n".join(case_entries))
    f.write("\n]\n")
    f.write(stage_block)
    f.write(pin_block)
    f.write(annexb_block)

print(f"wrote {out_path}: {len(cases)} cases, full stage dump = {full_dump_case!r}"
      + (", + Annex B derived case" if annexb_block else ""))
