# MBS-GO

Geometrical Optics (GO) light scattering by faceted particles. Computes Mueller matrices for ice crystals and other dielectric particles using beam-splitting ray tracing with Jones matrices and Fresnel coefficients.

## Build

Requires: `g++`, `make`, `qt5-qmake` (or `qmake-qt5`), `libqt5gui5`

```bash
cd pro
qmake MBS-GO.pro   # or qmake-qt5 MBS-GO.pro
make -j$(nproc)
```

Binary: `bin/mbs`

## Quick Start

```bash
# Hex column H=20 D=10 µm, ice, 10 reflections, backscattering
./bin/mbs -p 1 20 10 --ri 1.3116 0 -n 10 --fixed 0 0 \
          -w 0.532 --grid 0 180 1 1 --close

# Random orientation average (100x100 grid)
./bin/mbs -p 1 20 10 --ri 1.3116 0 -n 10 --random 100 100 \
          -w 0.532 --grid 0 180 1 1 --close
```

## Command-Line Reference

### Required

| Argument | Description |
|----------|-------------|
| `-p TYPE ...` | Particle type and dimensions (see below) |
| `--ri RE IM` | Complex refractive index (real and imaginary parts) |
| `-n N` | Number of internal reflections |
| `--grid ...` | Scattering angle grid (see below) |
| `--fixed β γ` or `--random Nβ Nγ` | Orientation: fixed (degrees) or random averaging |

### Optional

| Argument | Description |
|----------|-------------|
| `-w λ` | Wavelength in µm (required for `--adda` and `--abs`) |
| `-o NAME` | Output file/directory prefix (default: `M`) |
| `--close` | Exit after computation (no "press key" prompt) |
| `--abs` | Account for absorption (uses Im part of `--ri`) |
| `--incoh` | Incoherent scattering (no interference between beams) |
| `-r RATIO` | Beam area restriction ratio for splitting (default: 100) |
| `--log SEC` | Progress output interval in seconds |
| `--forced_convex` | Force convex particle tracing |
| `--forced_nonconvex` | Force non-convex particle tracing |
| `--tr FILE` | Compute only trajectories listed in file |
| `--all` | Compute all trajectories (with `--tr`) |

### Particle Types (`-p`)

| Type | Syntax | Description |
|------|--------|-------------|
| 1 | `-p 1 H D` | Hexagonal column (height H, diameter D in µm) |
| 2 | `-p 2 H D` | Bullet |
| 3 | `-p 3 H D [S]` | Bullet rosette |
| 4 | `-p 4 _ _ S` | Droxtal |
| 5 | `-p 5 L` | Cube (edge length L in µm) |
| 10 | `-p 10 H D C` | Concave hexagonal (cavity depth C) |
| 11 | `-p 11 H D` | Tilted hexagonal |
| 12 | `-p 12 H D N` | Hexagonal aggregate (N elements) |

### Scattering Grid (`--grid`)

3 parameters: `--grid R Nφ Nθ` — backscattering cone of radius R° from 180°

4 parameters: `--grid θ₁ θ₂ Nφ Nθ` — scattering angle range θ₁ to θ₂ (degrees)

Example: `--grid 0 180 1 1` — full range, forward to backward

## ADDA Bridge

Generate initial field approximation for [ADDA](https://github.com/adda-code/adda) (Discrete Dipole Approximation). MBS-GO computes the internal GO field and exports it as ADDA's `-init_field read` files.

### Usage

```bash
# Step 1: MBS-GO generates shape + field files
./bin/mbs -p 1 20 10 --ri 1.3116 0 -n 10 --fixed 0 0 \
          -w 0.532 --grid 0 180 1 1 --adda --dpl 10 --close -o test

# Step 2: ADDA uses them as initial field
adda -shape read test_shape.dat -dpl 10 -m 1.3116 0 \
     -prop 0 0 -1 -init_field read test_fieldY.dat test_fieldX.dat
```

### ADDA Mode Flags

| Argument | Description |
|----------|-------------|
| `--adda` | Enable ADDA mode (requires `--fixed` and `-w`) |
| `--dpl N` | Dipoles per wavelength (default: 10) |
| `--noinit` | Output only shape file, skip field files (for x₀=0 baseline) |
| `--norefl` | Direct refracted wave only, no internal reflections |
| `--rmax R` | Max \|Rs\|,\|Rp\| for analytical reflections (default: 0.25) |
| `--go` | Use GO-traced reflected beams instead of analytical |
| `--fp` | Fabry-Perot reflection model (anti-parallel facets only) |
| `--diffr` | Kirchhoff diffraction for GO reflections (with `--go`) |
| `--nf N` | Min Fresnel number for GO beams (with `--go`, default: 1.0) |
| `--rscale α` | Amplitude scaling for GO reflections (with `--go`, default: 1.0) |

### Output Files

| File | Format | Content |
|------|--------|---------|
| `{name}_shape.dat` | DDSCAT6 | Dipole coordinates |
| `{name}_fieldY.dat` | ADDA init_field | Y-polarization field (x y z \|E\|² Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i) |
| `{name}_fieldX.dat` | ADDA init_field | X-polarization field |

### Reflection Modes

**Default (analytical reflections)** — best overall. For each illuminated entry facet, traces a refracted ray to the opposite exit facet, computes Fresnel reflection Rs, Rp, and adds a reflected plane wave. Filters: anti-parallel facets only (nn_dot < -0.5), |R| < 0.25, no TIR. Helps at normal incidence (RE₀ reduced 8-22%), never hurts at oblique.

**`--go`** — uses actual GO-traced beam segments (nActs=1). Accurate amplitude/phase from Jones matrices but sharp beam boundaries cause artifacts at oblique incidence.

**`--fp`** — legacy Fabry-Perot model. Only works for anti-parallel facet pairs with |R| < 0.35.

**`--norefl`** — only the direct refracted wave, no reflections at all.

### Important Notes

- ADDA `-eps` is an **exponent**: `-eps 5` means ε = 10⁻⁵. Using `-eps 1e-5` gives ε ≈ 1 (no iterations!)
- Particle is rotated, not light: `-prop 0 0 -1` is always correct for ADDA
- ADDA convention: incPolY = (0,1,0), incPolX = (-1,0,0) when prop = (0,0,-1)

See [MBS-GO/ADDA_BRIDGE.md](MBS-GO/ADDA_BRIDGE.md) for detailed technical documentation, test results, and physical conventions.

## License

GNU General Public License v3.0
