# MBS-GO

Geometrical Optics (GO) light scattering by faceted particles. Computes Mueller matrices for ice crystals and other dielectric particles using beam-splitting ray tracing with Jones matrices and Fresnel coefficients.

## Build

Requires: `g++`, `make`

```bash
make -j$(nproc)
```

Binary: `bin/mbs`
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

### Flags

| Argument | Description |
|----------|-------------|
| `-p TYPE ...` | Particle type and dimensions (see below) |
| `--ri RE IM` | Complex refractive index (real and imaginary parts) |
| `-n N` | Number of internal reflections |
| `-dpl N` | Dipoles per lambda (default: 10)  |
| `--fixed β γ` | Particle orientation: zenith and azimuth angles (degrees) |
| `-w λ` | Wavelength in µm (required for `--adda`) |
| --adda |  Internal field mode calculation for ADDA

### Optional

| Argument | Description |
|----------|-------------|
| `--norefl` | Direct refracted wave only, no reflections |
| `--fp` | Legacy Fabry-Perot reflection model (anti-parallel facets only) |
| `--jmax J` | Max Jones matrix norm filter (default: no filter) |
| `--diffr` | Kirchhoff diffraction weighting for reflected beams |
| `--nf N` | Min Fresnel number for reflected beams (default: 0.0) |
| `--rscale α` | Amplitude scaling for reflected beams (default: 1.0) |
| `--goi` | Incoherent reflected beam accumulation |
| `--noinit` | Output only shape file, skip field files (for x₀=0 baseline) |
| `-r RATIO` | Beam area restriction ratio for splitting (default: 100) |
| `--forced_convex` | Force convex particle tracing |
| `--forced_nonconvex` | Force non-convex particle tracing |
| `--tr FILE` | Compute only trajectories listed in file |
| `--all` | Compute all trajectories (with `--tr`) |
| `-o NAME` | Output file/directory prefix (default: `M`) |

### Particle Types (`-p`)

| Type | Syntax | Description |
|------|--------|-------------|
| 1 | `-p 1 H D` | Hexagonal column (height H, diameter D in µm) |
| 2 | `-p 2 H D` | Bullet |
| 3 | `-p 3 H D [S]` | Bullet rosette |
| 4 | `-p 4 alpha1 alpha2 R` | Droxtal |
| 5 | `-p 5 L` | Cube (edge length L in µm) |
| 6 | `-p 6 D nLat nLon` | UV-sphere (diameter D, latitude × longitude tessellation) |
| 10 | `-p 10 H D C` | Concave hexagonal (cavity depth C) |
| 11 | `-p 11 H D` | Tilted hexagonal |
| 12 | `-p 12 H D N` | Hexagonal aggregate (N elements) |

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
- Particle is rotated, not light: `-prop 0 0 -1` is always correct for ADDA
- ADDA convention: incPolY = (0,1,0), incPolX = (-1,0,0) when prop = (0,0,-1)

See [MBS-GO/ADDA_BRIDGE_MANUAL.md](MBS-GO/ADDA_BRIDGE_MANUAL.md) for detailed technical documentation, test results, and physical conventions.

## License

GNU General Public License v3.0
