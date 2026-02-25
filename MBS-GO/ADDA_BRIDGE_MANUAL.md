# MBS-GO → ADDA Bridge: Initial Field for DDA Solver

## 1. Overview

The MBS-GO → ADDA bridge computes an initial approximation of the internal electromagnetic field inside a crystal particle using Geometrical Optics (GO), and passes it to the ADDA DDA solver via `-init_field read`. This reduces the initial residual norm (RE_000).

The quality metric is **RE_000** — the initial residual norm reported by ADDA. RE_000 < 1.0 means the GO-based initial field is better than zero (the default). Lower is better.

**Method: Per-facet Plane Wave + GO-traced reflected beams** — achieves RE_000 = 0.33–0.43 across all incidence angles for hexagonal columns, reducing the initial residual by 57–67%.

# 2. Usage

## 2.1 Command Line

```bash
# Generate ADDA input files
./bin/mbs -p <type> <height> <diameter> \
          --ri <Re_n> <Im_n> \
          -n <reflections> \
          --grid 0 180 1 1 \
          --fixed <beta> <gamma> \
          -w <wavelength_um> \
          --adda --dpl <dpl>

# Run ADDA with the generated initial field
adda -shape read M_shape.dat \
     -dpl <dpl> -m <Re_n> <Im_n> \
     -prop 0 0 -1 \
     -init_field read M_fieldY.dat M_fieldX.dat
```

## 2.2 Key Flags

| Flag | Description |
|------|-------------|
| `--adda` | Enable ADDA bridge mode |
| `--dpl N` | Dipoles per lambda (grid resolution, default: 10) |
| `--fixed beta gamma` | Incidence angles in degrees |
| `-n N` | Number of GO internal reflections |

## 2.3 Reflection Control Flags

By default, GO-traced reflected beams are accumulated with all segments up to `nActs = N` (set by `-n`) and no amplitude filtering.

| Flag | Description |
|------|-------------|
| `--norefl` | Skip all reflections |
| `--fp` | Use legacy Fabry-Perot model instead of GO beams |
| `--jmax J` | Max Jones matrix norm filter — skip beams with \|J\| > J (default: no filter) |
| `--nf N` | Min Fresnel number for GO beams (default: 0.0 = no filter) |
| `--rscale α` | Amplitude scaling factor for reflected beams (default: 1.0) |
| `--diffr` | Enable Kirchhoff diffraction weighting for reflected beams |
| `--goi` | Incoherent accumulation of reflected beams (intensity sum, not field sum) |

## 2.4 Other Flags

| Flag | Description |
|------|-------------|
| `--noinit` | Output only shape file, skip field files (for x₀=0 baseline) |

## 2.5 Output Files

- `M_shape.dat` — ADDA geometry file (dipole positions)
- `M_fieldY.dat` — Initial E-field for Y-polarization
- `M_fieldX.dat` — Initial E-field for X-polarization

## 2.6 Example

```bash
# Hexagonal column 10x5 um, n=1.3116, beta=30 deg, dpl=10
./bin/mbs -p 1 10 5 --ri 1.3116 0.0 -n 5 \
          --grid 0 180 1 1 --fixed 30 0 \
          -w 0.532 --adda --dpl 10

# Run ADDA
adda -shape read M_shape.dat -dpl 10 \
     -m 1.3116 0 -prop 0 0 -1 \
     -init_field read M_fieldY.dat M_fieldX.dat
```

# 3. Per-facet PW Method

## 3.1 Algorithm

The method assigns each internal dipole to one illuminated entry facet and fills it with a refracted plane wave:

**Step 1.** For each illuminated facet f (where cosθ_i = k_inc · n_in > 0):

- Compute refracted direction via Snell's law: k_refr = η·k_inc + C·n, where η = 1/n, C = cosθ_t − η·cosθ_i

- Compute complex Snell/Fresnel: γ = sqrt(m² − sin²θ_i), then Ts = 2cosθ_i / (cosθ_i + γ), Tp = 2m·cosθ_i / (m²cosθ_i + γ)

- Decompose incident polarization into s/p components and compute refracted E-field amplitude (complex for absorbing m)

**Step 2.** For each dipole at position r:

- Back-project along −k_refr to the facet plane
- Check if projection falls inside the facet polygon
- If yes: assign to this facet, compute field with complex phase: E(r) = A_f · exp[ik(k_inc · Δr_inc + (γ − cosθ_i)(n · Δr))], where Δr = r − r_ref. For absorbing m, Im(γ) > 0 provides exponential decay along the facet normal.
- If no facet matches: fallback to the primary facet (largest cosθ_i)

## 3.2 Key Properties

- **No overlap**: each dipole is assigned to exactly one facet via back-projection
- **Correct Fresnel**: s and p transmission coefficients applied separately
- **Size-independent**: RE_000 does not depend on `dpl` (confirmed for dpl=10 and dpl=20)
- **Phase-consistent**: uses k_inc · r for external path, n·k_refr · Δr for internal path

# 4. GO-traced Reflected Beams (Default)

## 4.1 Method

By default, the bridge uses GO-traced beam segments captured during ray tracing. These segments carry full Jones matrices from complex Fresnel reflections/transmissions along the beam path. For absorbing particles (Im(m) > 0), the Fresnel coefficients use complex Snell's law (cosθ_t = sqrt(1 − m²sin²θ_i)) and absorption is applied as exp(−k·κ·d) along internal paths.

For each reflected beam segment:

1. **Polygon cross-section**: the beam's exit polygon on the last facet defines its spatial extent
2. **Jones matrix**: encodes amplitude and phase from all Fresnel interactions along the path
3. **For each dipole**: check if it projects inside the beam polygon, compute E-field from Jones matrix with proper phase

## 4.2 Default Parameters

| Parameter | Default | Flag to change |
|-----------|---------|----------------|
| Max reflections (nActs) | = `-n` value | (set by `-n`) |
| Jones norm filter | no filter (1e10) | `--jmax J` |
| Fresnel number filter | no filter (0.0) | `--nf N` |
| Amplitude scaling | 1.0 | `--rscale α` |
| Diffraction weighting | off | `--diffr` |
| Accumulation mode | coherent | `--goi` for incoherent |

The philosophy is **minimal filtering by default** — all GO-traced beams are included as-is. Use the filter flags to restrict if needed.

## 4.3 Legacy: Fabry-Perot (`--fp`)

The `--fp` flag activates a simplified analytical reflection model:

- For each illuminated entry facet, trace refracted ray to the nearest exit facet
- Requires entry and exit normals to be anti-parallel (n_entry · n_exit < −0.95)
- Skips TIR cases and large reflections (|Rs| > 0.25 or |Rp| > 0.25)
- Applies Fabry-Perot multi-bounce factor: F = 1/(1 − R²_avg · exp(2iknd))

This mode is conservative and only helps for nearly anti-parallel facet pairs (e.g., basal faces of hexagonal columns at normal incidence).

# 5. Performance Results

## 5.1 RE_000 vs Incidence Angle (Hexagonal Column)

n = 1.3116, λ = 0.532 µm, dpl=10.

| Angle β | 10x5 µm (Y/X) | 20x10 µm (Y/X) |
|:-:|:-:|:-:|
| 0° | 0.427 / 0.426 | 0.357 / 0.356 |
| 10° | 0.383 / 0.377 | --- |
| 20° | 0.364 / 0.343 | --- |
| 30° | 0.374 / 0.336 | 0.317 / 0.274 |
| 45° | 0.375 / 0.328 | --- |
| 60° | 0.402 / 0.358 | 0.304 / 0.311 |

RE_000 improves with particle size: at 20x10 µm (~40λ) it reaches 0.27–0.36 (64–73% better than zero). GO approximation becomes more accurate as the size parameter grows.

## 5.2 Comparison of Methods (Hex Column, 30°)

| Method | RE_000 Y | RE_000 X | Notes |
|--------|:--------:|:--------:|-------|
| Zero field (default) | 1.000 | 1.000 | Baseline |
| Single PW (one direction) | 0.912 | 0.907 | Worst non-trivial |
| **Per-facet PW + GO refl** | **0.374** | **0.336** | **Best** |
| GO nActs=0 (exit-polygon) | 0.634 | 0.588 | 25% overlap |
| GO all beams | 0.789 | 0.801 | Overlap dominates |

## 5.3 Different Particle Shapes

All particles: n = 1.3116, λ = 0.532 µm, dpl=10.

| Particle | Angle | RE_000 Y/X | Notes |
|----------|:-----:|:----------:|-------|
| Hex column (type 1) | 0° | 0.427 / 0.426 | best at 0° |
| | 30° | 0.374 / 0.336 | |
| | 60° | 0.402 / 0.358 | |
| Concave hex (type 10) | 0° | 0.481 / 0.479 | cavity adds complexity |
| | 30° | 0.404 / 0.362 | |
| | 60° | 0.406 / 0.363 | |
| Droxtal (type 4) | 0° | 0.411 / 0.405 | 20 facets, quasi-sphere |
| | 30° | 0.349 / 0.333 | |
| | 60° | 0.314 / 0.316 | best shape overall |
| Bullet (type 2) | 0° | 1.374 / 1.385 | worse than zero! |
| | 30° | 0.782 / 0.751 | pyramidal top |
| | 60° | 0.398 / 0.334 | good, side facets |

**Key finding**: the method requires large, well-illuminated facets. Bullet crystals at normal incidence have a pyramidal top with small, steeply tilted facets --- per-facet PW fails (RE_000 > 1).

## 5.4 Size Scaling

RE_000 depends on physical particle size but NOT on grid resolution (dpl):

| Particle | dpl | Dipoles | RE_000 Y/X (30°) |
|----------|:---:|:-------:|:-:|
| 10x5 µm | 10 | 1.08M | 0.374 / 0.336 |
| 10x5 µm | 20 | 8.64M | 0.375 / 0.336 |
| 20x10 µm | 10 | 8.64M | 0.317 / 0.274 |

Doubling the particle size (at fixed dpl) improves RE_000 by ~0.06. Doubling dpl (at fixed size) does not change RE_000. The improvement with size confirms that GO becomes more accurate at larger size parameters.

# 6. Approaches Tested and Rejected

## 6.1 UTD/Kirchhoff Edge Diffraction for Entry Facets

Attempted to smooth entry facet boundaries using Kirchhoff polygon aperture weighting (`KirchhoffPolygonWeight`). Two variants tested:

- **Multi-facet coherent sum**: each dipole receives weighted contributions from ALL facets. Result: destructive interference, RE_000 = 2.1 (much worse than baseline).
- **Single-facet smoothing**: each dipole assigned to one facet, but with Fresnel edge smoothing instead of binary in/out. Result: RE_000 = 2.1 again — Fresnel scale sqrt(λ_eff · z_prop) was comparable to facet size, corrupting plane-wave phase structure.

**Conclusion**: Kirchhoff edge diffraction is unsuitable for facets with size < ~10λ.

## 6.2 Alternative Phase Models

| Flag | Method | RE_000 Y/X (hex 30°) | Status |
|------|--------|:-----:|--------|
| (default) | Complex Snell/Fresnel per-facet | 0.329 / 0.317 | **Best** |
| `--wkb` | Phase along incDir inside medium | 0.333 / 0.331 | Removed |
| `--pwkb` | Pure WKB (no refraction) | 0.454 / 0.551 | Removed |
| `--sph` | Spherical wavefront from facet center | worse | Removed |
| `--coh` | Coherent sum from all facets | worse | Removed |

## 6.3 Iteration Count

GO initial field does **not** reduce ADDA iteration count for QMR solver (the default). The system matrix condition number depends on particle shape and refractive index, not on the initial field. BiCGStab solver shows ~3.5% fewer iterations at m >= 1.5, but not for ice (n = 1.3116).

# 7. Code Details

## 7.1 Source Files

| File | Role |
|------|------|
| `src/adda/ADDAField.h` | Data structures: `DipoleField`, `InternalBeamSegment` |
| `src/adda/ADDAField.cpp` | Core implementation |
| `src/main.cpp` | ADDA mode entry point |
| `src/scattering/ScatteringConvex.cpp` | GO segment capture (convex) |
| `src/scattering/ScatteringNonConvex.cpp` | GO segment capture (non-convex) |

## 7.2 Key Functions

- `BuildDipoleGrid()` --- generates cubic grid, tests each point for particle membership
- `FillUncoveredPerFacet(incidentDir)` --- per-facet PW fill with complex Fresnel/Snell
- `AccumulateReflectedBeams(segments, incDir, maxJN, diffr, minNF, rScale, maxActs, incoh)` --- GO-traced reflected beam accumulation (default)
- `AddPerFacetReflection(incidentDir)` --- legacy Fabry-Perot reflection (`--fp`)
- `WriteGeometryFile()` / `WriteFieldFileY()` / `WriteFieldFileX()` --- ADDA output
- `DiagnoseGOvsPW(incidentDir)` --- diagnostic: compares GO field vs independent PW

## 7.3 Supported Particle Types

| Type | Shape | Convex | ADDA tested |
|:----:|-------|:------:|:-----------:|
| 1 | Hexagonal column | yes | yes |
| 2 | Bullet | yes | yes |
| 3 | Bullet rosette | no | --- |
| 4 | Droxtal | yes | yes |
| 5 | Cube | yes | yes |
| 6 | UV-sphere | yes | yes |
| 10 | Concave hexagonal | no | yes |
| 12 | Hexagonal aggregate | no | --- |

## 7.4 Coordinate System

- Particle is **rotated** by (β, γ, 0) before GO tracing
- ADDA uses fixed propagation direction `prop = (0, 0, -1)`
- All field coordinates are in the rotated particle frame

## 7.5 Phase Convention

φ = k · [k_inc · Δr_inc + (γ − cosθ_i)(n · Δr)]

where γ = sqrt(m² − sin²θ_i) (complex for absorbing m), n is the inward facet normal, and Δr = r − r_ref. For transparent m, this reduces to k_inc · r_proj + n · k_refr · (r − r_proj). For absorbing m, Im(γ) provides natural amplitude decay along the facet normal direction.

# 8. Known Issues and Bugs Fixed

1. **Snell's law sign bug** (fixed): was C = η·cosθ_i − cosθ_t, corrected to C = cosθ_t − η·cosθ_i

2. **E_par sign bug** (fixed): `e_par = cross(d, e_perp)` gives correct right-handed basis

3. **X-pol sign bug** (fixed): X-pol input is (0, −1) in (e_perp, e_par) basis, requiring sign flip

4. **Jones matrix s/p swap** (fixed): `AccumulateBeamContribution` had input columns swapped — Y-pol (pure s-input) was using Column 1 (p-input) instead of Column 2. Fix: Y-pol uses `(J.m22, J.m12)` = (s→s, s→p), X-pol uses `(-J.m21, -J.m11)` = (-(p→s), -(p→p)). At normal incidence this caused r_p (negative) to be applied instead of r_s (positive), making reflected beams subtract instead of add.

5. **GO beam overlap** (fundamental limitation): exit-polygon containment doesn't partition the interior --- multiple segments cover the same dipoles. Mitigated by using per-facet PW for the direct refracted wave and GO beams only for reflections.

6. **FAR_ZONE_DISTANCE** (fixed): was set to 10000, adding a large constant to optical path that corrupted phase for large wavelengths. Set to 0, then removed entirely.

# 9. Future Directions

- Test at larger size parameters (dpl >= 50) where GO approximation improves
- Implement Voronoi-like beam partitioning to eliminate reflected beam overlap
- Investigate boundary smoothing methods that preserve plane-wave phase structure
