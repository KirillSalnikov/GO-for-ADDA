# MBS-GO

Geometrical Optics (GO) light scattering by faceted particles. Computes Mueller matrices for ice crystals and other dielectric particles using beam-splitting ray tracing with Jones matrices and Fresnel coefficients.

## Build

Requires: `g++`, `make`

```bash
make -j$(nproc)
```

Binary: `bin/mbs`

## Quick Start

```bash
# Hex column H=20 D=10 µm, ice, 10 reflections
./bin/mbs -p 1 20 10 --ri 1.3116 0 -n 10 --fixed 0 0 \
          -w 0.532 --grid 0 180 1 1
# UV-sphere D=50 µm, absorbing particle (m = 1.5 + 0.01i)
./bin/mbs -p 6 50 20 20 --ri 1.5 0.01 -n 10 --fixed 0 0 \
          -w 0.532 --grid 0 180 1 1```

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
| `-w λ` | Wavelength in µm (required for `--adda`) |
| `-o NAME` | Output file/directory prefix (default: `M`) |
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
| 6 | `-p 6 D nLat nLon` | UV-sphere (diameter D, latitude × longitude tessellation) |
| 10 | `-p 10 H D C` | Concave hexagonal (cavity depth C) |
| 11 | `-p 11 H D` | Tilted hexagonal |
| 12 | `-p 12 H D N` | Hexagonal aggregate (N elements) |

### Refractive Index (`--ri`)

`--ri RE IM` — complex refractive index m = RE + i·IM. Both parts are always used: RE controls refraction (Snell's law), IM controls absorption (exponential decay along internal paths). Fresnel coefficients use full complex Snell's law.

Examples: `--ri 1.3116 0` (ice, transparent), `--ri 1.5 0.01` (absorbing glass), `--ri 1.3 0.1` (strongly absorbing).

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
          -w 0.532 --grid 0 180 1 1 --adda --dpl 10 -o test

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
| `--norefl` | Direct refracted wave only, no reflections |
| `--fp` | Legacy Fabry-Perot reflection model (anti-parallel facets only) |
| `--jmax J` | Max Jones matrix norm filter (default: no filter) |
| `--diffr` | Kirchhoff diffraction weighting for reflected beams |
| `--nf N` | Min Fresnel number for reflected beams (default: 0.0) |
| `--rscale α` | Amplitude scaling for reflected beams (default: 1.0) |
| `--goi` | Incoherent reflected beam accumulation |

### Output Files

| File | Format | Content |
|------|--------|---------|
| `{name}_shape.dat` | DDSCAT6 | Dipole coordinates |
| `{name}_fieldY.dat` | ADDA init_field | Y-polarization field (x y z \|E\|² Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i) |
| `{name}_fieldX.dat` | ADDA init_field | X-polarization field |

### Reflection Modes

**Default (GO-traced beams)** — uses actual GO-traced beam segments with full Jones matrices and complex Fresnel coefficients. All reflected beams up to `nActs = -n` are included with no amplitude filtering. Use `--jmax`, `--nf` to restrict if needed. Absorption (Im(m) > 0) is handled automatically via complex Snell's law and exponential decay along internal paths.

**`--fp`** — legacy Fabry-Perot analytical model. Only works for anti-parallel facet pairs with small Fresnel coefficients.

**`--norefl`** — only the direct refracted wave, no reflections at all.

### Important Notes

- ADDA `-eps` is an **exponent**: `-eps 5` means ε = 10⁻⁵. Using `-eps 1e-5` gives ε ≈ 1 (no iterations!)
- Particle is rotated, not light: `-prop 0 0 -1` is always correct for ADDA
- ADDA convention: incPolY = (0,1,0), incPolX = (-1,0,0) when prop = (0,0,-1)

See [MBS-GO/ADDA_BRIDGE_MANUAL.md](MBS-GO/ADDA_BRIDGE_MANUAL.md) for detailed technical documentation, test results, and physical conventions.

## License

GNU General Public License v3.0
