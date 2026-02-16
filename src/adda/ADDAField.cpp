#include "ADDAField.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

/// Compute Fresnel integrals C(x) = ∫₀ˣ cos(πt²/2)dt, S(x) = ∫₀ˣ sin(πt²/2)dt
static void FresnelCS(double x, double &outC, double &outS)
{
    double ax = fabs(x);
    if (ax < 1.6)
    {
        // Power series: C+iS = Σ_{n=0} (iπ/2)^n x^{2n+1} / (n!(2n+1))
        double x2 = ax * ax;
        double re = 0, im = 0;
        double tr = ax, ti = 0;  // term: starts real = x
        double piH = M_PI / 2.0;
        for (int n = 0; n <= 20; ++n)
        {
            re += tr / (2*n + 1);
            im += ti / (2*n + 1);
            // Multiply term by (iπ/2) * x² / (n+1)
            double f = piH * x2 / (n + 1);
            double tr2 = -ti * f;  // real part of (tr+i*ti) * i*f
            double ti2 =  tr * f;  // imag part
            tr = tr2;
            ti = ti2;
        }
        outC = re;
        outS = im;
    }
    else
    {
        // Asymptotic: C = 0.5 + f*sin - g*cos, S = 0.5 - f*cos - g*sin
        // Rational approximation for f(x), g(x) (Abramowitz & Stegun §7.3)
        double f = (1.0 + 0.926*ax) / (2.0 + 1.792*ax + 3.104*ax*ax);
        double g = 1.0 / (2.0 + 4.142*ax + 3.492*ax*ax + 6.670*ax*ax*ax);
        double pix2h = M_PI * ax * ax / 2.0;
        double sn = sin(pix2h), cs = cos(pix2h);
        outC = 0.5 + f*sn - g*cs;
        outS = 0.5 - f*cs - g*sn;
    }
    if (x < 0) { outC = -outC; outS = -outS; }
}

/// Fresnel half-plane diffraction weight (complex)
/// d = signed distance to edge (>0 illuminated, <0 shadow)
/// z = propagation distance from aperture
/// lambda_eff = wavelength in the medium
static complex FresnelEdgeWeight(double d, double z, double lambda_eff)
{
    if (z < 1e-10) return complex(d > 0 ? 1.0 : 0.0, 0.0);
    double w = d * sqrt(2.0 / (lambda_eff * z));
    double C, S;
    FresnelCS(w, C, S);
    // U/U_0 = 0.5 + (1+i)/2 * (C - iS) = 0.5 + (C+S)/2 + i(C-S)/2
    return complex(0.5 + (C + S) / 2.0, (C - S) / 2.0);
}

/// Kirchhoff diffraction weight from a polygon aperture.
/// Uses Green's theorem to reduce the 2D Fresnel integral to 1D edge integrals:
///   ∫∫_poly exp(iπ(X²+Y²)/2) dXdY = ∮ Q(X,Y) dY
/// where Q = exp(iπY²/2) · [(1+i)/2 + C(X) + iS(X)]
/// Each edge is integrated with 8-point Gauss-Legendre quadrature.
/// Returns U/U₀ = (-i/2) · ∮ Q dY  (normalized to 1 for infinite aperture).
static complex KirchhoffPolygonWeight(
    const double *poly_u, const double *poly_v, int nVert,
    double obs_u, double obs_v,
    double z, double lambda_eff)
{
    if (z < 1e-10 || nVert < 3)
        return complex(1.0, 0.0);

    double scale = sqrt(2.0 / (lambda_eff * z));

    // 8-point Gauss-Legendre on [0,1]
    static const int NG = 8;
    static const double gp[NG] = {
        0.01985507175123, 0.10166676129319, 0.23723379504184, 0.40828267875218,
        0.59171732124782, 0.76276620495816, 0.89833323870681, 0.98014492824877
    };
    static const double gw[NG] = {
        0.05061426814519, 0.11119051722669, 0.15685332293895, 0.18134189168918,
        0.18134189168918, 0.15685332293895, 0.11119051722669, 0.05061426814519
    };

    complex result(0.0, 0.0);

    for (int ei = 0; ei < nVert; ++ei)
    {
        int ej = (ei + 1) % nVert;

        // Fresnel-scaled coordinates relative to observation point
        double X1 = (poly_u[ei] - obs_u) * scale;
        double Y1 = (poly_v[ei] - obs_v) * scale;
        double X2 = (poly_u[ej] - obs_u) * scale;
        double Y2 = (poly_v[ej] - obs_v) * scale;

        double dY = Y2 - Y1;
        if (fabs(dY) < 1e-15) continue;  // dY=0 → ∫Q dY = 0 (correct by Green's thm)
        double dX = X2 - X1;

        complex edgeSum(0.0, 0.0);
        for (int g = 0; g < NG; ++g)
        {
            double t = gp[g];
            double X = X1 + t * dX;
            double Y = Y1 + t * dY;

            double C, S;
            FresnelCS(X, C, S);

            double piY2h = M_PI * Y * Y / 2.0;
            double cosY = cos(piY2h), sinY = sin(piY2h);

            // exp(iπY²/2) · [(1+i)/2 + C(X) + iS(X)]
            double ar = 0.5 + C;
            double ai = 0.5 + S;
            double re = cosY * ar - sinY * ai;
            double im = sinY * ar + cosY * ai;

            edgeSum += complex(re, im) * gw[g];
        }

        result += edgeSum * dY;
    }

    // U/U₀ = (1/(2i)) · result = (-i/2) · result
    return complex(0.0, -0.5) * result;
}

ADDAFieldComputer::ADDAFieldComputer(Particle *particle, double wavelength,
                                     const complex &ri, int dpl)
    : m_particle(particle),
      m_wavelength(wavelength),
      m_ri(ri),
      m_segCount(0),
      m_zMax(0),
      m_boxX(0), m_boxY(0), m_boxZ(0)
{
    m_gridSpacing = wavelength / dpl;
    m_k = 2.0 * M_PI / wavelength;
}

bool ADDAFieldComputer::IsInsideParticle(const Point3f &p) const
{
    if (!m_particle->isConcave)
    {
        // For convex particle: point inside iff on inner side of every facet
        for (int i = 0; i < m_particle->nFacets; ++i)
        {
            const Point3f &n = m_particle->facets[i].in_normal;
            double dp = (double)p.cx * n.cx + (double)p.cy * n.cy
                      + (double)p.cz * n.cz + n.d_param;
            if (dp < -1e-6)
                return false;
        }
        return true;
    }

    // For non-convex particle: ray casting along +Z axis
    // Count intersections of ray from p in +Z direction with facet triangles
    int crossings = 0;
    Point3f rayDir(0, 0, 1);

    for (int i = 0; i < m_particle->nFacets; ++i)
    {
        const Facet &f = m_particle->facets[i];
        const Point3f &n = f.ex_normal;

        // Ray-plane intersection: t = -(n·p + d) / (n·rayDir)
        double nDotRay = (double)n.cz;  // n · (0,0,1)
        if (fabs(nDotRay) < 1e-12)
            continue;  // ray parallel to facet

        double nDotP = (double)p.cx * n.cx + (double)p.cy * n.cy
                     + (double)p.cz * n.cz + n.d_param;
        double t = -nDotP / nDotRay;
        if (t < 1e-6)
            continue;  // intersection behind ray origin

        // Intersection point
        double hitZ = p.cz + t;
        Point3f hit((float)p.cx, (float)p.cy, (float)hitZ);

        // Check if hit point is inside the facet polygon
        // Triangulate facet as fan from vertex 0
        for (int j = 1; j < f.nVertices - 1; ++j)
        {
            const Point3f &v0 = f.arr[0];
            const Point3f &v1 = f.arr[j];
            const Point3f &v2 = f.arr[j + 1];

            // Barycentric test in XY plane (since ray is along Z)
            double dx1 = v1.cx - v0.cx, dy1 = v1.cy - v0.cy;
            double dx2 = v2.cx - v0.cx, dy2 = v2.cy - v0.cy;
            double dxP = hit.cx - v0.cx, dyP = hit.cy - v0.cy;

            double det = dx1 * dy2 - dx2 * dy1;
            if (fabs(det) < 1e-12)
                continue;

            double u = (dxP * dy2 - dx2 * dyP) / det;
            double v = (dx1 * dyP - dxP * dy1) / det;

            if (u >= -1e-6 && v >= -1e-6 && (u + v) <= 1.0 + 1e-6)
            {
                ++crossings;
                break;  // count each facet at most once
            }
        }
    }

    return (crossings % 2) == 1;
}

void ADDAFieldComputer::BuildDipoleGrid()
{
    // Compute bounding box of particle
    double minX = 1e30, maxX = -1e30;
    double minY = 1e30, maxY = -1e30;
    double minZ = 1e30, maxZ = -1e30;

    for (int i = 0; i < m_particle->nFacets; ++i)
    {
        for (int j = 0; j < m_particle->facets[i].nVertices; ++j)
        {
            const Point3f &v = m_particle->facets[i].arr[j];
            if (v.cx < minX) minX = v.cx;
            if (v.cx > maxX) maxX = v.cx;
            if (v.cy < minY) minY = v.cy;
            if (v.cy > maxY) maxY = v.cy;
            if (v.cz < minZ) minZ = v.cz;
            if (v.cz > maxZ) maxZ = v.cz;
        }
    }

    m_zMax = maxZ;

    // Grid dimensions (add margin of 1 cell on each side)
    m_boxX = (int)ceil((maxX - minX) / m_gridSpacing) + 2;
    m_boxY = (int)ceil((maxY - minY) / m_gridSpacing) + 2;
    m_boxZ = (int)ceil((maxZ - minZ) / m_gridSpacing) + 2;

    // Grid origin: center grid on particle center
    double centerX = (minX + maxX) / 2.0;
    double centerY = (minY + maxY) / 2.0;
    double centerZ = (minZ + maxZ) / 2.0;

    double x0 = centerX - (m_boxX - 1) * m_gridSpacing / 2.0;
    double y0 = centerY - (m_boxY - 1) * m_gridSpacing / 2.0;
    double z0 = centerZ - (m_boxZ - 1) * m_gridSpacing / 2.0;

    m_dipoles.clear();

    for (int iz = 0; iz < m_boxZ; ++iz)
    {
        for (int iy = 0; iy < m_boxY; ++iy)
        {
            for (int ix = 0; ix < m_boxX; ++ix)
            {
                Point3f p((float)(x0 + ix * m_gridSpacing),
                          (float)(y0 + iy * m_gridSpacing),
                          (float)(z0 + iz * m_gridSpacing));

                if (IsInsideParticle(p))
                {
                    DipoleField dip;
                    dip.x = p.cx;
                    dip.y = p.cy;
                    dip.z = p.cz;
                    dip.ix = ix;
                    dip.iy = iy;
                    dip.iz = iz;
                    m_dipoles.push_back(dip);
                }
            }
        }
    }

    cout << "ADDA grid: " << m_boxX << "x" << m_boxY << "x" << m_boxZ
         << ", dipoles inside: " << m_dipoles.size()
         << ", spacing: " << m_gridSpacing << " um" << endl;
}

bool ADDAFieldComputer::IsInsidePolygon(const Point3f &proj,
                                         const Point3f *polyArr, int nVert,
                                         const Point3f &normal) const
{
    // Convex polygon containment: point inside if on same side of all edges
    for (int i = 0; i < nVert; ++i)
    {
        int j = (i + 1) % nVert;
        Point3f edge = polyArr[j] - polyArr[i];
        Point3f toPoint = proj - polyArr[i];
        Point3f cross = CrossProduct(edge, toPoint);
        double dp = (double)cross.cx * normal.cx
                  + (double)cross.cy * normal.cy
                  + (double)cross.cz * normal.cz;
        if (dp < -1e-4)
            return false;
    }
    return true;
}

double ADDAFieldComputer::MinEdgeDistance(const Point3f &proj,
                                          const Point3f *polyArr, int nVert) const
{
    double minDist = 1e30;
    for (int i = 0; i < nVert; ++i)
    {
        int j = (i + 1) % nVert;
        // Edge vector and length
        double ex = polyArr[j].cx - polyArr[i].cx;
        double ey = polyArr[j].cy - polyArr[i].cy;
        double ez = polyArr[j].cz - polyArr[i].cz;
        double edgeLen2 = ex*ex + ey*ey + ez*ez;
        if (edgeLen2 < 1e-20) continue;

        // Project point onto edge line: t = dot(proj - v[i], edge) / |edge|²
        double dx = proj.cx - polyArr[i].cx;
        double dy = proj.cy - polyArr[i].cy;
        double dz = proj.cz - polyArr[i].cz;
        double t = (dx*ex + dy*ey + dz*ez) / edgeLen2;
        t = (t < 0.0) ? 0.0 : (t > 1.0 ? 1.0 : t);

        // Distance from proj to closest point on edge segment
        double px = dx - t*ex;
        double py = dy - t*ey;
        double pz = dz - t*ez;
        double dist = sqrt(px*px + py*py + pz*pz);
        if (dist < minDist) minDist = dist;
    }
    return minDist;
}

void ADDAFieldComputer::AccumulateBeamContribution(const InternalBeamSegment &seg)
{
    // Beam local polarization basis
    Point3d e_perp(seg.polarizationBasis.cx,
                   seg.polarizationBasis.cy,
                   seg.polarizationBasis.cz);
    Point3d d(seg.direction.cx, seg.direction.cy, seg.direction.cz);
    Point3d e_par = CrossProductD(d, e_perp);
    double len = LengthD(e_par);
    if (len < 1e-12) return;
    e_par = e_par / len;

    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    // Compute polygon normal from exit polygon vertices for containment test
    Point3f polyNormal(0, 0, 0);
    if (seg.polyNVertices >= 3)
    {
        Point3f e1 = seg.polyArr[1] - seg.polyArr[0];
        Point3f e2 = seg.polyArr[2] - seg.polyArr[0];
        polyNormal = CrossProduct(e1, e2);
        float pnLen = (float)Length(polyNormal);
        if (pnLen > 1e-10f)
            polyNormal = polyNormal / (double)pnLen;
    }

    for (size_t di = 0; di < m_dipoles.size(); ++di)
    {
        DipoleField &dip = m_dipoles[di];

        // First-write-wins: skip dipoles already covered by another beam
        if (dip.nContributions > 0)
            continue;

        Point3f r((float)dip.x, (float)dip.y, (float)dip.z);

        // Distance from beam entry plane to dipole along beam direction
        double t = (double)seg.direction.cx * r.cx
                 + (double)seg.direction.cy * r.cy
                 + (double)seg.direction.cz * r.cz
                 + seg.front;

        // Check within segment bounds
        if (t < -1e-6 || t > seg.segmentLength + 1e-6)
            continue;

        // Project dipole to exit facet plane along beam direction
        Point3f exitNorm = seg.exitNormal;
        double dp_nr = (double)r.cx * exitNorm.cx
                     + (double)r.cy * exitNorm.cy
                     + (double)r.cz * exitNorm.cz
                     + seg.exitD;
        double dp_nd = (double)seg.direction.cx * exitNorm.cx
                     + (double)seg.direction.cy * exitNorm.cy
                     + (double)seg.direction.cz * exitNorm.cz;

        if (fabs(dp_nd) < 1e-10)
            continue;

        double t_exit = -dp_nr / dp_nd;
        Point3f r_proj(r.cx + seg.direction.cx * (float)t_exit,
                       r.cy + seg.direction.cy * (float)t_exit,
                       r.cz + seg.direction.cz * (float)t_exit);

        // Check if projection is inside the intersection polygon
        if (!IsInsidePolygon(r_proj, seg.polyArr, seg.polyNVertices, polyNormal))
            continue;

        // Compute phase from GO optical path
        double delta_OP = t * re_n;
        double total_OP = seg.opticalPath + delta_OP;
        complex phase = exp_im(m_k * total_OP);

        // Absorption
        if (im_n > 1e-15)
        {
            phase = phase * exp(-m_k * im_n * t);
        }

        ++dip.nContributions;

        // Jones matrix: first-write (=), not accumulate (+=)
        // e_perp = s-direction, e_par = p-direction
        // Y-polarization: pure s-input → (0, 1) in (p, s) basis
        complex Ep_Y = seg.J.m22 * phase;
        complex Epar_Y = seg.J.m12 * phase;
        dip.Ex_Y = Ep_Y * e_perp.x + Epar_Y * e_par.x;
        dip.Ey_Y = Ep_Y * e_perp.y + Epar_Y * e_par.y;
        dip.Ez_Y = Ep_Y * e_perp.z + Epar_Y * e_par.z;

        // X-polarization: pure -p-input → (-1, 0) in (p, s) basis
        complex Ep_X = -seg.J.m21 * phase;
        complex Epar_X = -seg.J.m11 * phase;
        dip.Ex_X = Ep_X * e_perp.x + Epar_X * e_par.x;
        dip.Ey_X = Ep_X * e_perp.y + Epar_X * e_par.y;
        dip.Ez_X = Ep_X * e_perp.z + Epar_X * e_par.z;
    }
}

void ADDAFieldComputer::WriteGeometryFile(const string &filename) const
{
    ofstream file(filename);
    if (!file.is_open())
    {
        cerr << "ERROR: Cannot open " << filename << " for writing" << endl;
        return;
    }

    file << "MBS-GO generated geometry; grid_spacing=" << m_gridSpacing
         << " um; wavelength=" << m_wavelength << " um" << endl;
    file << m_dipoles.size() << " = NAT" << endl;
    file << "1 0 0 = A_1 vector" << endl;
    file << "0 1 0 = A_2 vector" << endl;
    file << "1.0 1.0 1.0 = lattice spacings (d_x,d_y,d_z)/d" << endl;
    file << "JA  IX  IY  IZ ICOMP1 ICOMP2 ICOMP3" << endl;

    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        const DipoleField &dip = m_dipoles[i];
        file << (i + 1) << " "
             << dip.ix << " " << dip.iy << " " << dip.iz << " "
             << "1 1 1" << endl;
    }

    file.close();
    cout << "Wrote geometry: " << filename << " (" << m_dipoles.size() << " dipoles)" << endl;
}

void ADDAFieldComputer::WriteFieldFileY(const string &filename) const
{
    ofstream file(filename);
    if (!file.is_open())
    {
        cerr << "ERROR: Cannot open " << filename << " for writing" << endl;
        return;
    }

    file << scientific << setprecision(10);
    file << "x y z |E|^2 Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i" << endl;

    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        const DipoleField &dip = m_dipoles[i];
        double E2 = norm(dip.Ex_Y) + norm(dip.Ey_Y) + norm(dip.Ez_Y);

        file << dip.x << " " << dip.y << " " << dip.z << " "
             << E2 << " "
             << real(dip.Ex_Y) << " " << imag(dip.Ex_Y) << " "
             << real(dip.Ey_Y) << " " << imag(dip.Ey_Y) << " "
             << real(dip.Ez_Y) << " " << imag(dip.Ez_Y) << endl;
    }

    file.close();
    cout << "Wrote Y-pol field: " << filename << endl;
}

void ADDAFieldComputer::WriteFieldFileX(const string &filename) const
{
    ofstream file(filename);
    if (!file.is_open())
    {
        cerr << "ERROR: Cannot open " << filename << " for writing" << endl;
        return;
    }

    file << scientific << setprecision(10);
    file << "x y z |E|^2 Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i" << endl;

    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        const DipoleField &dip = m_dipoles[i];
        double E2 = norm(dip.Ex_X) + norm(dip.Ey_X) + norm(dip.Ez_X);

        file << dip.x << " " << dip.y << " " << dip.z << " "
             << E2 << " "
             << real(dip.Ex_X) << " " << imag(dip.Ex_X) << " "
             << real(dip.Ey_X) << " " << imag(dip.Ey_X) << " "
             << real(dip.Ez_X) << " " << imag(dip.Ez_X) << endl;
    }

    file.close();
    cout << "Wrote X-pol field: " << filename << endl;
}

void ADDAFieldComputer::FillUncoveredWithPlaneWave(const Point3f &refractedDir,
                                                     const Point3f &incidentDir,
                                                     double Ts, double Tp)
{
    // Fill uncovered dipoles (nContributions==0) with refracted plane wave.
    // Uses the same phase and polarization logic as ApplyUniformBeam.
    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    double dx = (double)refractedDir.cx;
    double dy = (double)refractedDir.cy;
    double dz = (double)refractedDir.cz;
    double ix = (double)incidentDir.cx;
    double iy = (double)incidentDir.cy;
    double iz = (double)incidentDir.cz;

    // s/p polarization basis from plane of incidence
    double sx = iy * dz - iz * dy;
    double sy = iz * dx - ix * dz;
    double sz = ix * dy - iy * dx;
    double slen = sqrt(sx*sx + sy*sy + sz*sz);

    double EY_x, EY_y, EY_z;
    double EX_x, EX_y, EX_z;

    if (slen < 1e-6)
    {
        double T = 0.5 * (Ts + Tp);
        EY_x = 0; EY_y = T; EY_z = 0;
        EX_x = -T; EX_y = 0; EX_z = 0;
    }
    else
    {
        sx /= slen; sy /= slen; sz /= slen;
        double pix = sy * iz - sz * iy;
        double piy = sz * ix - sx * iz;
        double piz = sx * iy - sy * ix;
        double pilen = sqrt(pix*pix + piy*piy + piz*piz);
        pix /= pilen; piy /= pilen; piz /= pilen;
        double prx = sy * dz - sz * dy;
        double pry = sz * dx - sx * dz;
        double prz = sx * dy - sy * dx;
        double prlen = sqrt(prx*prx + pry*pry + prz*prz);
        prx /= prlen; pry /= prlen; prz /= prlen;

        double Eys = sy, Eyp = piy;
        EY_x = Ts * Eys * sx + Tp * Eyp * prx;
        EY_y = Ts * Eys * sy + Tp * Eyp * pry;
        EY_z = Ts * Eys * sz + Tp * Eyp * prz;

        double Exs = -sx, Exp = -pix;
        EX_x = Ts * Exs * sx + Tp * Exp * prx;
        EX_y = Ts * Exs * sy + Tp * Exp * pry;
        EX_z = Ts * Exs * sz + Tp * Exp * prz;
    }

    static const double FAR_ZONE_DISTANCE = 10000.0;
    double externalOP = FAR_ZONE_DISTANCE + iz * m_zMax;
    double dRefDotRef = dz * m_zMax;

    int uncovered = 0;
    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        DipoleField &dip = m_dipoles[i];
        if (dip.nContributions == 0)
        {
            double t = dx * dip.x + dy * dip.y + dz * dip.z - dRefDotRef;
            double total_OP = externalOP + re_n * t;
            complex phase = exp_im(m_k * total_OP);

            if (im_n > 1e-15)
            {
                double tAbs = (t > 0) ? t : 0;
                phase = phase * exp(-m_k * im_n * tAbs);
            }

            dip.Ex_Y = EY_x * phase;
            dip.Ey_Y = EY_y * phase;
            dip.Ez_Y = EY_z * phase;
            dip.Ex_X = EX_x * phase;
            dip.Ey_X = EX_y * phase;
            dip.Ez_X = EX_z * phase;
            ++uncovered;
        }
    }
    cout << "Dipoles covered by GO: " << (m_dipoles.size() - uncovered)
         << ", fallback to refracted plane wave: " << uncovered << endl;
}

void ADDAFieldComputer::FillUncoveredPerFacet(const Point3f &incidentDir)
{
    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    double ix = (double)incidentDir.cx;
    double iy = (double)incidentDir.cy;
    double iz = (double)incidentDir.cz;

    // Collect illuminated facets with their refraction parameters
    struct FacetRefr {
        int id;
        double cosI;
        Point3f refDir;
        double Ts, Tp;
        // E-field amplitude vectors for Y-pol and X-pol
        double EY_x, EY_y, EY_z;
        double EX_x, EX_y, EX_z;
        // Phase reference: FAR_ZONE + incDir . facetCenter
        double phaseRef;
    };

    static const double FAR_ZONE_DISTANCE = 10000.0;
    vector<FacetRefr> facetList;

    for (int f = 0; f < m_particle->nFacets; ++f)
    {
        const Facet &fac = m_particle->facets[f];
        double nx = (double)fac.in_normal.cx;
        double ny = (double)fac.in_normal.cy;
        double nz = (double)fac.in_normal.cz;
        double cosI = ix*nx + iy*ny + iz*nz;
        if (cosI < 1e-6)
            continue;

        // Snell's law: refracted direction
        double sinIsq = 1.0 - cosI*cosI;
        double cosT = sqrt(fmax(0.0, 1.0 - sinIsq/(re_n*re_n)));
        double eta = 1.0 / re_n;
        double C = cosT - eta * cosI;
        double rdx = eta * ix + C * nx;
        double rdy = eta * iy + C * ny;
        double rdz = eta * iz + C * nz;
        double rlen = sqrt(rdx*rdx + rdy*rdy + rdz*rdz);
        rdx /= rlen; rdy /= rlen; rdz /= rlen;

        // Fresnel transmission coefficients
        double Ts = 2.0 * cosI / (cosI + re_n * cosT);
        double Tp = 2.0 * cosI / (re_n * cosI + cosT);

        // s/p polarization basis from plane of incidence (incDir, refDir)
        double sx = iy*rdz - iz*rdy;
        double sy = iz*rdx - ix*rdz;
        double sz = ix*rdy - iy*rdx;
        double slen = sqrt(sx*sx + sy*sy + sz*sz);

        double EY_x, EY_y, EY_z;
        double EX_x, EX_y, EX_z;

        if (slen < 1e-6)
        {
            // Normal incidence: Ts ~ Tp
            double T = 0.5 * (Ts + Tp);
            EY_x = 0; EY_y = T; EY_z = 0;
            EX_x = -T; EX_y = 0; EX_z = 0;
        }
        else
        {
            sx /= slen; sy /= slen; sz /= slen;
            // p_inc = cross(s, incDir)
            double pix_ = sy*iz - sz*iy;
            double piy_ = sz*ix - sx*iz;
            double piz_ = sx*iy - sy*ix;
            double pilen = sqrt(pix_*pix_ + piy_*piy_ + piz_*piz_);
            pix_ /= pilen; piy_ /= pilen; piz_ /= pilen;
            // p_refr = cross(s, refDir)
            double prx = sy*rdz - sz*rdy;
            double pry = sz*rdx - sx*rdz;
            double prz = sx*rdy - sy*rdx;
            double prlen = sqrt(prx*prx + pry*pry + prz*prz);
            prx /= prlen; pry /= prlen; prz /= prlen;

            // Y-pol: E_inc = (0,1,0)
            double Eys = sy;
            double Eyp = piy_;
            EY_x = Ts*Eys*sx + Tp*Eyp*prx;
            EY_y = Ts*Eys*sy + Tp*Eyp*pry;
            EY_z = Ts*Eys*sz + Tp*Eyp*prz;

            // X-pol: E_inc = (-1,0,0)
            double Exs = -sx;
            double Exp = -pix_;
            EX_x = Ts*Exs*sx + Tp*Exp*prx;
            EX_y = Ts*Exs*sy + Tp*Exp*pry;
            EX_z = Ts*Exs*sz + Tp*Exp*prz;
        }

        FacetRefr fr;
        fr.id = f;
        fr.cosI = cosI;
        fr.refDir = Point3f((float)rdx, (float)rdy, (float)rdz);
        fr.Ts = Ts;
        fr.Tp = Tp;
        fr.EY_x = EY_x; fr.EY_y = EY_y; fr.EY_z = EY_z;
        fr.EX_x = EX_x; fr.EX_y = EX_y; fr.EX_z = EX_z;
        // Phase ref: external OP to facet center
        double cx_ = (double)fac.center.cx;
        double cy_ = (double)fac.center.cy;
        double cz_ = (double)fac.center.cz;
        fr.phaseRef = FAR_ZONE_DISTANCE + ix*cx_ + iy*cy_ + iz*cz_;
        facetList.push_back(fr);
    }

    // Sort by cosI descending (primary facet first for fallback)
    for (size_t i = 0; i < facetList.size(); ++i)
        for (size_t j = i+1; j < facetList.size(); ++j)
            if (facetList[j].cosI > facetList[i].cosI)
                swap(facetList[i], facetList[j]);

    // If all facets are steeply tilted (e.g., pyramidal top at normal incidence),
    // use incident direction for phase propagation instead of refracted direction.
    // Tilted refracted PWs from individual facets create wrong phase patterns
    // deep inside the crystal; a straight-through PW is a better approximation.
    bool pyramidMode = (!facetList.empty() && facetList[0].cosI < 0.6);

    cout << "Per-facet refraction: " << facetList.size() << " illuminated facets";
    if (!facetList.empty())
        cout << " (primary cosI=" << facetList[0].cosI << ")";
    if (pyramidMode)
        cout << " [pyramid mode: using incDir for phase]";
    cout << endl;

    int uncovered = 0;
    int perFacetFilled = 0;
    int fallbackFilled = 0;

    if (pyramidMode)
    {
        // Pyramid mode: all facets are steeply tilted, per-facet PW creates
        // wrong phase patterns. Use a single PW along incDir with
        // area-weighted average Fresnel amplitude.
        double totalWeight = 0;
        double avgEY_x = 0, avgEY_y = 0, avgEY_z = 0;
        double avgEX_x = 0, avgEX_y = 0, avgEX_z = 0;
        double avgCx = 0, avgCy = 0, avgCz = 0;
        for (size_t fi = 0; fi < facetList.size(); ++fi) {
            const FacetRefr &fr = facetList[fi];
            const Facet &fac = m_particle->facets[fr.id];
            double w = fr.cosI; // projected area weight
            totalWeight += w;
            avgEY_x += w * fr.EY_x; avgEY_y += w * fr.EY_y; avgEY_z += w * fr.EY_z;
            avgEX_x += w * fr.EX_x; avgEX_y += w * fr.EX_y; avgEX_z += w * fr.EX_z;
            avgCx += w * (double)fac.center.cx;
            avgCy += w * (double)fac.center.cy;
            avgCz += w * (double)fac.center.cz;
        }
        avgEY_x /= totalWeight; avgEY_y /= totalWeight; avgEY_z /= totalWeight;
        avgEX_x /= totalWeight; avgEX_y /= totalWeight; avgEX_z /= totalWeight;
        avgCx /= totalWeight; avgCy /= totalWeight; avgCz /= totalWeight;
        double avgPhaseRef = FAR_ZONE_DISTANCE + ix*avgCx + iy*avgCy + iz*avgCz;

        cout << "  Pyramid mode: avg E_Y=(" << avgEY_x << "," << avgEY_y << ","
             << avgEY_z << "), ref=(" << avgCx << "," << avgCy << "," << avgCz << ")" << endl;

        for (size_t i = 0; i < m_dipoles.size(); ++i)
        {
            DipoleField &dip = m_dipoles[i];
            if (dip.nContributions > 0)
                continue;

            ++uncovered;
            double t_int = ix*(dip.x - avgCx) + iy*(dip.y - avgCy) + iz*(dip.z - avgCz);
            double total_OP = avgPhaseRef + re_n * t_int;
            complex phase = exp_im(m_k * total_OP);

            if (im_n > 1e-15)
            {
                double tAbs = (t_int > 0) ? t_int : 0;
                phase = phase * exp(-m_k * im_n * tAbs);
            }

            dip.Ex_Y = avgEY_x * phase;
            dip.Ey_Y = avgEY_y * phase;
            dip.Ez_Y = avgEY_z * phase;
            dip.Ex_X = avgEX_x * phase;
            dip.Ey_X = avgEX_y * phase;
            dip.Ez_X = avgEX_z * phase;
            dip.assignedFacet = facetList[0].id;
            ++perFacetFilled;
        }
    }
    else
    {
    // Normal mode: per-facet assignment via back-projection
    // Always set assignedFacet (for FP reflection); fill field only for uncovered dipoles
    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        DipoleField &dip = m_dipoles[i];
        bool goCovered = (dip.nContributions > 0);
        if (!goCovered)
            ++uncovered;

        bool found = false;

        for (size_t fi = 0; fi < facetList.size(); ++fi)
        {
            const FacetRefr &fr = facetList[fi];
            const Facet &fac = m_particle->facets[fr.id];

            // Back-project dipole along -incidentDir to facet plane
            double dp_nr = (double)dip.x * fac.in_normal.cx
                         + (double)dip.y * fac.in_normal.cy
                         + (double)dip.z * fac.in_normal.cz
                         + fac.in_normal.d_param;
            double dp_nd = ix * fac.in_normal.cx
                         + iy * fac.in_normal.cy
                         + iz * fac.in_normal.cz;
            if (fabs(dp_nd) < 1e-10)
                continue;

            double t_back = -dp_nr / dp_nd;
            if (t_back > 1e-6)
                continue;

            Point3f proj((float)(dip.x + ix * t_back),
                         (float)(dip.y + iy * t_back),
                         (float)(dip.z + iz * t_back));

            Point3f polyNormal(0, 0, 0);
            if (fac.nVertices >= 3)
            {
                Point3f e1 = fac.arr[1] - fac.arr[0];
                Point3f e2 = fac.arr[2] - fac.arr[0];
                polyNormal = CrossProduct(e1, e2);
                float pnLen = (float)Length(polyNormal);
                if (pnLen > 1e-10f)
                    polyNormal = polyNormal / (double)pnLen;
            }

            if (!IsInsidePolygon(proj, fac.arr, fac.nVertices, polyNormal))
                continue;

            // Always set assignedFacet (needed for FP reflection)
            dip.assignedFacet = fr.id;
            found = true;

            if (!goCovered)
            {
                // Fill field with per-facet PW using exact entry point
                double rdx = (double)fr.refDir.cx;
                double rdy = (double)fr.refDir.cy;
                double rdz = (double)fr.refDir.cz;
                double px = (double)proj.cx;
                double py = (double)proj.cy;
                double pz = (double)proj.cz;
                double extOP = FAR_ZONE_DISTANCE + ix*px + iy*py + iz*pz;
                double t_int = rdx*(dip.x - px) + rdy*(dip.y - py) + rdz*(dip.z - pz);
                double total_OP = extOP + re_n * t_int;
                complex phase = exp_im(m_k * total_OP);

                if (im_n > 1e-15)
                {
                    double tAbs = (t_int > 0) ? t_int : 0;
                    phase = phase * exp(-m_k * im_n * tAbs);
                }

                dip.Ex_Y = fr.EY_x * phase;
                dip.Ey_Y = fr.EY_y * phase;
                dip.Ez_Y = fr.EY_z * phase;
                dip.Ex_X = fr.EX_x * phase;
                dip.Ey_X = fr.EX_y * phase;
                dip.Ez_X = fr.EX_z * phase;
                ++perFacetFilled;
            }
            break;
        }

        if (!found && !facetList.empty())
        {
            // Fallback: use primary facet (largest cosI)
            const FacetRefr &fr = facetList[0];
            dip.assignedFacet = fr.id;

            if (!goCovered)
            {
                const Facet &fac = m_particle->facets[fr.id];
                double rdx = (double)fr.refDir.cx;
                double rdy = (double)fr.refDir.cy;
                double rdz = (double)fr.refDir.cz;
                double fb_nr = (double)dip.x * fac.in_normal.cx
                             + (double)dip.y * fac.in_normal.cy
                             + (double)dip.z * fac.in_normal.cz
                             + fac.in_normal.d_param;
                double fb_nd = ix * fac.in_normal.cx
                             + iy * fac.in_normal.cy
                             + iz * fac.in_normal.cz;
                double fb_t = (fabs(fb_nd) > 1e-10) ? -fb_nr / fb_nd : 0.0;
                double px = dip.x + ix * fb_t;
                double py = dip.y + iy * fb_t;
                double pz = dip.z + iz * fb_t;
                double extOP = FAR_ZONE_DISTANCE + ix*px + iy*py + iz*pz;
                double t_int = rdx*(dip.x - px) + rdy*(dip.y - py) + rdz*(dip.z - pz);
                double total_OP = extOP + re_n * t_int;
                complex phase = exp_im(m_k * total_OP);

                if (im_n > 1e-15)
                {
                    double tAbs = (t_int > 0) ? t_int : 0;
                    phase = phase * exp(-m_k * im_n * tAbs);
                }

                dip.Ex_Y = fr.EY_x * phase;
                dip.Ey_Y = fr.EY_y * phase;
                dip.Ez_Y = fr.EY_z * phase;
                dip.Ex_X = fr.EX_x * phase;
                dip.Ey_X = fr.EX_y * phase;
                dip.Ez_X = fr.EX_z * phase;
                ++fallbackFilled;
            }
        }
    }
    } // end normal mode

    cout << "Dipoles covered by GO: " << (m_dipoles.size() - uncovered)
         << ", per-facet plane wave: " << perFacetFilled
         << ", primary-facet fallback: " << fallbackFilled << endl;
}

void ADDAFieldComputer::AddPerFacetReflection(const Point3f &incidentDir)
{
    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    double ix = (double)incidentDir.cx;
    double iy = (double)incidentDir.cy;
    double iz = (double)incidentDir.cz;

    static const double FAR_ZONE_DISTANCE = 10000.0;

    // For each illuminated entry facet, compute reflected beam parameters
    struct FacetRefl {
        int entryId;
        int exitId;           // exit facet that refracted beam hits
        double cosI_entry;    // cos(incidence) at entry
        Point3f refDir;       // refracted direction (into crystal)
        Point3f reflDir;      // reflected direction (after bounce off exit facet)
        double Rs, Rp;        // Fresnel reflection coefficients at exit (amplitude)
        double Ts_entry, Tp_entry;  // Fresnel transmission at entry
        double phaseRef;      // phase reference for reflected wave
        double exitPt_x, exitPt_y, exitPt_z;  // actual exit point (not facet center!)
        complex fpFactor;     // Fabry-Perot multi-bounce: 1/(1 - R²*exp(2iknd))
        // Reflected E-field amplitude vectors
        double EY_x, EY_y, EY_z;
        double EX_x, EX_y, EX_z;
    };

    vector<FacetRefl> reflList;

    for (int f = 0; f < m_particle->nFacets; ++f)
    {
        const Facet &entry = m_particle->facets[f];
        double nx = (double)entry.in_normal.cx;
        double ny = (double)entry.in_normal.cy;
        double nz = (double)entry.in_normal.cz;
        double cosI = ix*nx + iy*ny + iz*nz;
        if (cosI < 1e-6)
            continue;

        // Entry refraction (Snell's law)
        double sinIsq = 1.0 - cosI*cosI;
        double cosT = sqrt(fmax(0.0, 1.0 - sinIsq/(re_n*re_n)));
        double eta = 1.0 / re_n;
        double C = cosT - eta * cosI;
        double rdx = eta*ix + C*nx;
        double rdy = eta*iy + C*ny;
        double rdz = eta*iz + C*nz;
        double rlen = sqrt(rdx*rdx + rdy*rdy + rdz*rdz);
        rdx /= rlen; rdy /= rlen; rdz /= rlen;

        // Entry Fresnel transmission
        double Ts_entry = 2.0*cosI / (cosI + re_n*cosT);
        double Tp_entry = 2.0*cosI / (re_n*cosI + cosT);

        // Find exit facet: trace refracted beam from entry center
        int exitId = -1;
        double minT = 1e30;
        for (int g = 0; g < m_particle->nFacets; ++g)
        {
            if (g == f) continue;
            const Facet &exitFac = m_particle->facets[g];
            // Ray-plane intersection: (refDir . exitNormal) * t = -(r0 . exitNormal + d)
            double dnr = rdx*exitFac.in_normal.cx + rdy*exitFac.in_normal.cy + rdz*exitFac.in_normal.cz;
            if (dnr > -1e-6) continue;  // exit facet: refracted beam hits outer side (in_normal points inward, beam hits from inside)
            double dp = entry.center.cx*exitFac.in_normal.cx
                      + entry.center.cy*exitFac.in_normal.cy
                      + entry.center.cz*exitFac.in_normal.cz
                      + exitFac.in_normal.d_param;
            double t = -dp / dnr;
            if (t > 1e-6 && t < minT)
            {
                minT = t;
                exitId = g;
            }
        }
        if (exitId < 0) continue;

        const Facet &exitFac = m_particle->facets[exitId];

        // Only use reflection for truly anti-parallel facet pairs (basal faces)
        double nn_dot = nx*exitFac.in_normal.cx + ny*exitFac.in_normal.cy + nz*exitFac.in_normal.cz;
        if (nn_dot > -0.95) continue;  // not anti-parallel → skip

        // Internal incidence angle at exit facet
        double enx = (double)exitFac.in_normal.cx;
        double eny = (double)exitFac.in_normal.cy;
        double enz = (double)exitFac.in_normal.cz;
        // cosI_exit = -refDir . in_normal (beam hits from inside, in_normal points inward)
        double cosI_exit = -(rdx*enx + rdy*eny + rdz*enz);
        if (cosI_exit < 1e-6) continue;

        // Check for total internal reflection
        double sinI_exit_sq = 1.0 - cosI_exit*cosI_exit;
        double sinT_exit_sq = sinI_exit_sq * re_n * re_n;  // n_crystal * sinI = n_air * sinT
        if (sinT_exit_sq > 1.0)
        {
            // Total internal reflection — full reflection
            // Rs = Rp = 1 (magnitude), but complex phase
            // For simplicity, use Fresnel formulas with complex cosT
            // Skip TIR for now — these beams are already captured by GO
            continue;
        }

        double cosT_exit = sqrt(1.0 - sinT_exit_sq);

        // Fresnel reflection at exit (crystal→air): amplitude coefficients
        // r_s = (n*cosI - cosT) / (n*cosI + cosT)
        // r_p = (cosI - n*cosT) / (cosI + n*cosT)
        double Rs = (re_n*cosI_exit - cosT_exit) / (re_n*cosI_exit + cosT_exit);
        double Rp = (cosI_exit - re_n*cosT_exit) / (cosI_exit + re_n*cosT_exit);

        // Skip large reflections — these worsen the approximation
        if (fabs(Rs) > 0.35 || fabs(Rp) > 0.35) continue;

        // Reflected direction at exit facet (specular reflection)
        // reflDir = refDir + 2*cosI_exit * in_normal (reflect off inner surface)
        // Actually: reflDir = refDir - 2*(refDir.exitNormal_out)*exitNormal_out
        // Since in_normal points inward and refDir . in_normal = -cosI_exit:
        // reflDir = refDir + 2*cosI_exit * in_normal
        double reflx = rdx + 2.0*cosI_exit*enx;
        double refly = rdy + 2.0*cosI_exit*eny;
        double reflz = rdz + 2.0*cosI_exit*enz;
        double refll = sqrt(reflx*reflx + refly*refly + reflz*reflz);
        reflx /= refll; refly /= refll; reflz /= refll;

        // s/p basis for decomposition
        // s_entry = cross(incDir, refDir) (same as in FillUncoveredPerFacet)
        double sx = iy*rdz - iz*rdy;
        double sy = iz*rdx - ix*rdz;
        double sz = ix*rdy - iy*rdx;
        double slen = sqrt(sx*sx + sy*sy + sz*sz);

        if (slen < 1e-6)
        {
            // Normal incidence: R_s ≈ -R_p at small angles
            // Use average |R|, with sign from R_s (physical convention)
            double R = Rs;  // at normal incidence, use s-convention (positive)
            double EY_amp = Ts_entry * R;
            double EX_amp = Ts_entry * R;

            FacetRefl fr;
            fr.entryId = f;
            fr.exitId = exitId;
            fr.cosI_entry = cosI;
            fr.refDir = Point3f((float)rdx, (float)rdy, (float)rdz);
            fr.reflDir = Point3f((float)reflx, (float)refly, (float)reflz);
            fr.Rs = Rs; fr.Rp = Rp;
            fr.Ts_entry = Ts_entry; fr.Tp_entry = Tp_entry;

            // Reflected Y-pol: same as direct but with R factor and reflected direction
            // At normal incidence, reflected E stays along Y
            fr.EY_x = 0; fr.EY_y = EY_amp; fr.EY_z = 0;
            fr.EX_x = -EX_amp; fr.EX_y = 0; fr.EX_z = 0;

            // Exit point = entry_center + minT * refDir (NOT exit facet center!)
            double cx_ = (double)entry.center.cx;
            double cy_ = (double)entry.center.cy;
            double cz_ = (double)entry.center.cz;
            fr.exitPt_x = cx_ + minT * rdx;
            fr.exitPt_y = cy_ + minT * rdy;
            fr.exitPt_z = cz_ + minT * rdz;
            fr.phaseRef = FAR_ZONE_DISTANCE + ix*cx_ + iy*cy_ + iz*cz_ + re_n * minT;

            // Fabry-Perot multi-bounce factor: 1/(1 - R²*exp(2iknd))
            double R2avg = 0.5*(Rs*Rs + Rp*Rp);
            complex rt = exp_im(2.0 * m_k * re_n * minT);
            if (im_n > 1e-15) rt = rt * exp(-2.0 * m_k * im_n * minT);
            double dr = 1.0 - R2avg * real(rt);
            double di = -R2avg * imag(rt);
            double dm2 = dr*dr + di*di;
            fr.fpFactor = complex(dr/dm2, -di/dm2);

            reflList.push_back(fr);
        }
        else
        {
            sx /= slen; sy /= slen; sz /= slen;
            // p_inc = cross(s_entry, incDir) — p for incident wave at entry
            double pix_ = sy*iz - sz*iy;
            double piy_ = sz*ix - sx*iz;
            double piz_ = sx*iy - sy*ix;
            double pilen = sqrt(pix_*pix_ + piy_*piy_ + piz_*piz_);
            pix_ /= pilen; piy_ /= pilen; piz_ /= pilen;
            // p_refracted = cross(s_entry, refDir) — p for refracted wave (entry basis)
            double pfx = sy*rdz - sz*rdy;
            double pfy = sz*rdx - sx*rdz;
            double pfz = sx*rdy - sy*rdx;
            double pflen = sqrt(pfx*pfx + pfy*pfy + pfz*pfz);
            pfx /= pflen; pfy /= pflen; pfz /= pflen;

            // Y-pol: E_inc = (0,1,0), decompose into entry s, p_inc
            double Eys = sy;   // dot((0,1,0), s_entry)
            double Eyp = piy_; // dot((0,1,0), p_inc)
            // Transmitted field inside crystal (3D vector):
            // E_inside = Ts*Eys * s_entry + Tp*Eyp * p_refracted
            double EinY_x = Ts_entry*Eys*sx + Tp_entry*Eyp*pfx;
            double EinY_y = Ts_entry*Eys*sy + Tp_entry*Eyp*pfy;
            double EinY_z = Ts_entry*Eys*sz + Tp_entry*Eyp*pfz;

            // X-pol: E_inc = (-1,0,0)
            double Exs = -sx;
            double Exp_ = -pix_;
            double EinX_x = Ts_entry*Exs*sx + Tp_entry*Exp_*pfx;
            double EinX_y = Ts_entry*Exs*sy + Tp_entry*Exp_*pfy;
            double EinX_z = Ts_entry*Exs*sz + Tp_entry*Exp_*pfz;

            // Exit s/p basis (plane of incidence at exit facet)
            // exitOutNormal = -exitInNormal
            double eox = -enx, eoy = -eny, eoz = -enz;
            double sex = rdy*eoz - rdz*eoy;  // s_exit = cross(refDir, exitOutNormal)
            double sey = rdz*eox - rdx*eoz;
            double sez = rdx*eoy - rdy*eox;
            double selen = sqrt(sex*sex + sey*sey + sez*sez);

            double rEY_x, rEY_y, rEY_z;
            double rEX_x, rEX_y, rEX_z;

            if (selen < 1e-6)
            {
                // Normal incidence at exit — Rs ≈ -Rp, use average R
                double R = (fabs(Rs) + fabs(Rp)) * 0.5;
                if (Rs > 0) R = R; else R = -R;
                rEY_x = R * EinY_x; rEY_y = R * EinY_y; rEY_z = R * EinY_z;
                rEX_x = R * EinX_x; rEX_y = R * EinX_y; rEX_z = R * EinX_z;
            }
            else
            {
                sex /= selen; sey /= selen; sez /= selen;
                // p for refracted beam at exit: p_exit_inc = cross(s_exit, refDir)
                double peix = sey*rdz - sez*rdy;
                double peiy = sez*rdx - sex*rdz;
                double peiz = sex*rdy - sey*rdx;
                double peilen = sqrt(peix*peix + peiy*peiy + peiz*peiz);
                peix /= peilen; peiy /= peilen; peiz /= peilen;

                // Decompose E_inside into exit s/p
                double EinY_se = EinY_x*sex + EinY_y*sey + EinY_z*sez;
                double EinY_pe = EinY_x*peix + EinY_y*peiy + EinY_z*peiz;
                double EinX_se = EinX_x*sex + EinX_y*sey + EinX_z*sez;
                double EinX_pe = EinX_x*peix + EinX_y*peiy + EinX_z*peiz;

                // Apply exit Fresnel reflection
                double ErY_se = Rs * EinY_se;
                double ErY_pe = Rp * EinY_pe;
                double ErX_se = Rs * EinX_se;
                double ErX_pe = Rp * EinX_pe;

                // p for reflected beam at exit: p_refl_exit = cross(s_exit, reflDir)
                double prex = sey*reflz - sez*refly;
                double prey = sez*reflx - sex*reflz;
                double prez = sex*refly - sey*reflx;
                double prelen = sqrt(prex*prex + prey*prey + prez*prez);
                prex /= prelen; prey /= prelen; prez /= prelen;

                // Reflected field (3D vector) = s_comp * s_exit + p_comp * p_refl_exit
                rEY_x = ErY_se*sex + ErY_pe*prex;
                rEY_y = ErY_se*sey + ErY_pe*prey;
                rEY_z = ErY_se*sez + ErY_pe*prez;
                rEX_x = ErX_se*sex + ErX_pe*prex;
                rEX_y = ErX_se*sey + ErX_pe*prey;
                rEX_z = ErX_se*sez + ErX_pe*prez;
            }

            FacetRefl fr;
            fr.entryId = f;
            fr.exitId = exitId;
            fr.cosI_entry = cosI;
            fr.refDir = Point3f((float)rdx, (float)rdy, (float)rdz);
            fr.reflDir = Point3f((float)reflx, (float)refly, (float)reflz);
            fr.Rs = Rs; fr.Rp = Rp;
            fr.Ts_entry = Ts_entry; fr.Tp_entry = Tp_entry;
            fr.EY_x = rEY_x; fr.EY_y = rEY_y; fr.EY_z = rEY_z;
            fr.EX_x = rEX_x; fr.EX_y = rEX_y; fr.EX_z = rEX_z;

            // Exit point = entry_center + minT * refDir
            double cx_ = (double)entry.center.cx;
            double cy_ = (double)entry.center.cy;
            double cz_ = (double)entry.center.cz;
            fr.exitPt_x = cx_ + minT * rdx;
            fr.exitPt_y = cy_ + minT * rdy;
            fr.exitPt_z = cz_ + minT * rdz;
            fr.phaseRef = FAR_ZONE_DISTANCE + ix*cx_ + iy*cy_ + iz*cz_ + re_n * minT;

            // Fabry-Perot multi-bounce factor: 1/(1 - R²*exp(2iknd))
            double R2avg = 0.5*(Rs*Rs + Rp*Rp);
            complex rt = exp_im(2.0 * m_k * re_n * minT);
            if (im_n > 1e-15) rt = rt * exp(-2.0 * m_k * im_n * minT);
            double dr = 1.0 - R2avg * real(rt);
            double di = -R2avg * imag(rt);
            double dm2 = dr*dr + di*di;
            fr.fpFactor = complex(dr/dm2, -di/dm2);

            reflList.push_back(fr);
        }
    }

    cout << "Per-facet reflection: " << reflList.size() << " facet pairs" << endl;

    // For each dipole, add reflected wave from its assigned entry facet
    int nAdded = 0;
    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        DipoleField &dip = m_dipoles[i];
        if (dip.assignedFacet < 0) continue;

        // Find the reflection entry for this facet
        const FacetRefl *fr = nullptr;
        for (size_t ri = 0; ri < reflList.size(); ++ri)
        {
            if (reflList[ri].entryId == dip.assignedFacet)
            {
                fr = &reflList[ri];
                break;
            }
        }
        if (!fr) continue;

        // Phase: phaseRef is at the actual exit point (entry_center + minT*refDir)
        // Reflected beam travels from exit point in reflDir
        // Distance from exit point to dipole along reflDir:
        double t_refl = (double)fr->reflDir.cx * (dip.x - fr->exitPt_x)
                       + (double)fr->reflDir.cy * (dip.y - fr->exitPt_y)
                       + (double)fr->reflDir.cz * (dip.z - fr->exitPt_z);

        // Dipole should be "behind" the exit facet in reflected direction
        // (reflected beam goes back into the crystal)
        if (t_refl < -1e-6) continue;

        double total_OP = fr->phaseRef + re_n * t_refl;
        complex phase = exp_im(m_k * total_OP);

        if (im_n > 1e-15)
        {
            double tAbs = (t_refl > 0) ? t_refl : 0;
            phase = phase * exp(-m_k * im_n * tAbs);
        }

        dip.Ex_Y += fr->EY_x * phase;
        dip.Ey_Y += fr->EY_y * phase;
        dip.Ez_Y += fr->EY_z * phase;
        dip.Ex_X += fr->EX_x * phase;
        dip.Ey_X += fr->EX_y * phase;
        dip.Ez_X += fr->EX_z * phase;

        // Apply Fabry-Perot multi-bounce factor to total field (direct + reflected)
        dip.Ex_Y *= fr->fpFactor;
        dip.Ey_Y *= fr->fpFactor;
        dip.Ez_Y *= fr->fpFactor;
        dip.Ex_X *= fr->fpFactor;
        dip.Ey_X *= fr->fpFactor;
        dip.Ez_X *= fr->fpFactor;
        ++nAdded;
    }

    cout << "Per-facet Fabry-Perot: " << reflList.size() << " pairs, " << nAdded << " dipoles" << endl;
}

void ADDAFieldComputer::AccumulateReflectedBeams(
    const std::vector<InternalBeamSegment> &segments,
    const Point3f &incidentDir, double maxJonesNorm)
{
    double re_n = real(m_ri);
    double im_n = imag(m_ri);
    static const double FAR_ZONE_DISTANCE = 10000.0;

    double iDirX = (double)incidentDir.cx;
    double iDirY = (double)incidentDir.cy;
    double iDirZ = (double)incidentDir.cz;

    int nProcessed = 0;
    int nSkippedNorm = 0;
    int nDipolesUpdated = 0;
    int nSmoothed = 0;
    int nEntryMissed = 0;

    for (size_t si = 0; si < segments.size(); ++si)
    {
        const InternalBeamSegment &seg = segments[si];
        if (seg.nActs != 1)
            continue;

        // Jones norm filter: skip beams with large |J| (e.g. near-TIR)
        double jNorm = sqrt(norm(seg.J.m11) + norm(seg.J.m12)
                          + norm(seg.J.m21) + norm(seg.J.m22));
        if (jNorm > maxJonesNorm)
        {
            ++nSkippedNorm;
            continue;
        }

        ++nProcessed;

        // Beam local polarization basis (same as AccumulateBeamContribution)
        Point3d e_perp(seg.polarizationBasis.cx,
                       seg.polarizationBasis.cy,
                       seg.polarizationBasis.cz);
        Point3d dDir(seg.direction.cx, seg.direction.cy, seg.direction.cz);
        Point3d e_par = CrossProductD(dDir, e_perp);
        double len = LengthD(e_par);
        if (len < 1e-12) continue;
        e_par = e_par / len;

        // Compute polygon normal for containment test
        Point3f polyNormal(0, 0, 0);
        if (seg.polyNVertices >= 3)
        {
            Point3f e1 = seg.polyArr[1] - seg.polyArr[0];
            Point3f e2 = seg.polyArr[2] - seg.polyArr[0];
            polyNormal = CrossProduct(e1, e2);
            float pnLen = (float)Length(polyNormal);
            if (pnLen > 1e-10f)
                polyNormal = polyNormal / (double)pnLen;
        }

        // Reflection facet normal (in_normal points into crystal)
        const Facet &reflFacet = m_particle->facets[seg.entryFacetId];
        double nBx = (double)reflFacet.in_normal.cx;
        double nBy = (double)reflFacet.in_normal.cy;
        double nBz = (double)reflFacet.in_normal.cz;
        double nBd = reflFacet.in_normal.d_param;

        // Recover d_inc (internal direction before reflection)
        // Reflection: d_refl = d_inc - 2*(d_inc·n_out)*n_out
        // With n_out = -in_normal: d_inc = d_refl - 2*(d_refl·in_normal)*in_normal
        double dDotN = dDir.x * nBx + dDir.y * nBy + dDir.z * nBz;
        double dIncX = dDir.x - 2.0 * dDotN * nBx;
        double dIncY = dDir.y - 2.0 * dDotN * nBy;
        double dIncZ = dDir.z - 2.0 * dDotN * nBz;

        double lambda_eff = m_wavelength / re_n;

        // Back-project exit polygon to reflection facet B for Kirchhoff integral.
        // 2D coordinates use e_perp/e_par (perpendicular to beam), NOT facet-plane axes.
        double polyB_u[MAX_VERTEX_NUM], polyB_v[MAX_VERTEX_NUM];
        double kb_ox = 0, kb_oy = 0, kb_oz = 0;
        bool kirchhoffReady = false;

        if (seg.polyNVertices >= 3 && fabs(dDotN) > 1e-10)
        {
            // Back-project first vertex to get origin on facet B
            double px0 = (double)seg.polyArr[0].cx;
            double py0 = (double)seg.polyArr[0].cy;
            double pz0 = (double)seg.polyArr[0].cz;
            double dp0 = px0*nBx + py0*nBy + pz0*nBz + nBd;
            double t0 = dp0 / dDotN;
            kb_ox = px0 - t0*dDir.x;
            kb_oy = py0 - t0*dDir.y;
            kb_oz = pz0 - t0*dDir.z;

            // Back-project all vertices and project onto beam-perpendicular axes
            for (int vi = 0; vi < seg.polyNVertices; ++vi)
            {
                double px = (double)seg.polyArr[vi].cx;
                double py = (double)seg.polyArr[vi].cy;
                double pz = (double)seg.polyArr[vi].cz;
                double dp_v = px*nBx + py*nBy + pz*nBz + nBd;
                double tv = dp_v / dDotN;
                double vx = px - tv*dDir.x - kb_ox;
                double vy = py - tv*dDir.y - kb_oy;
                double vz = pz - tv*dDir.z - kb_oz;
                polyB_u[vi] = vx*e_perp.x + vy*e_perp.y + vz*e_perp.z;
                polyB_v[vi] = vx*e_par.x  + vy*e_par.y  + vz*e_par.z;
            }

            // Ensure CCW orientation (Green's theorem requires it)
            double signedArea2D = 0;
            for (int vi = 0; vi < seg.polyNVertices; ++vi)
            {
                int vj = (vi + 1) % seg.polyNVertices;
                signedArea2D += polyB_u[vi]*polyB_v[vj] - polyB_u[vj]*polyB_v[vi];
            }
            if (signedArea2D < 0)
            {
                for (int vi = 0; vi < seg.polyNVertices / 2; ++vi)
                {
                    int vj = seg.polyNVertices - 1 - vi;
                    std::swap(polyB_u[vi], polyB_u[vj]);
                    std::swap(polyB_v[vi], polyB_v[vj]);
                }
            }
            kirchhoffReady = true;
        }

        int segUpdated = 0;
        for (size_t di = 0; di < m_dipoles.size(); ++di)
        {
            DipoleField &dip = m_dipoles[di];

            Point3f r((float)dip.x, (float)dip.y, (float)dip.z);

            // Check within segment bounds (beam direction distance)
            double t = (double)seg.direction.cx * r.cx
                     + (double)seg.direction.cy * r.cy
                     + (double)seg.direction.cz * r.cz
                     + seg.front;

            if (t < -1e-6 || t > seg.segmentLength + 1e-6)
                continue;

            // Check if dipole projects inside exit polygon
            Point3f exitNorm = seg.exitNormal;
            double dp_nr = (double)r.cx * exitNorm.cx
                         + (double)r.cy * exitNorm.cy
                         + (double)r.cz * exitNorm.cz
                         + seg.exitD;
            double dp_nd = (double)seg.direction.cx * exitNorm.cx
                         + (double)seg.direction.cy * exitNorm.cy
                         + (double)seg.direction.cz * exitNorm.cz;

            if (fabs(dp_nd) < 1e-10)
                continue;

            double t_exit = -dp_nr / dp_nd;
            Point3f r_proj(r.cx + seg.direction.cx * (float)t_exit,
                           r.cy + seg.direction.cy * (float)t_exit,
                           r.cz + seg.direction.cz * (float)t_exit);

            // Pre-filter uses 2D ray-casting (no normal dependency).
            // Quick back-project dipole to facet B → get 2D coords → check containment.
            if (kirchhoffReady && fabs(dDotN) > 1e-10)
            {
                double dp_rn = dip.x*nBx + dip.y*nBy + dip.z*nBz + nBd;
                double al = -dp_rn / dDotN;
                double rx = dip.x + al*dDir.x - kb_ox;
                double ry = dip.y + al*dDir.y - kb_oy;
                double rz = dip.z + al*dDir.z - kb_oz;
                double ou = rx*e_perp.x + ry*e_perp.y + rz*e_perp.z;
                double ov = rx*e_par.x  + ry*e_par.y  + rz*e_par.z;

                // 2D ray-casting point-in-polygon test
                int crossings = 0;
                for (int ei = 0; ei < seg.polyNVertices; ++ei)
                {
                    int ej = (ei + 1) % seg.polyNVertices;
                    double vi_v = polyB_v[ei], vj_v = polyB_v[ej];
                    if ((vi_v <= ov && vj_v > ov) || (vj_v <= ov && vi_v > ov))
                    {
                        double t_rc = (ov - vi_v) / (vj_v - vi_v);
                        if (ou < polyB_u[ei] + t_rc * (polyB_u[ej] - polyB_u[ei]))
                            ++crossings;
                    }
                }
                bool inside2D = (crossings & 1) != 0;

                if (!inside2D)
                {
                    // Outside: check minimum edge distance in 2D
                    double minEdge2D = 1e30;
                    for (int ei = 0; ei < seg.polyNVertices; ++ei)
                    {
                        int ej = (ei + 1) % seg.polyNVertices;
                        double eu = polyB_u[ej] - polyB_u[ei];
                        double ev = polyB_v[ej] - polyB_v[ei];
                        double el2 = eu*eu + ev*ev;
                        if (el2 < 1e-20) continue;
                        double du = ou - polyB_u[ei], dv = ov - polyB_v[ei];
                        double tc = (du*eu + dv*ev) / el2;
                        tc = (tc < 0) ? 0 : (tc > 1) ? 1 : tc;
                        double dx = du - tc*eu, dy = dv - tc*ev;
                        double d2 = dx*dx + dy*dy;
                        if (d2 < minEdge2D) minEdge2D = d2;
                    }
                    if (sqrt(minEdge2D) > 5.0 * lambda_eff)
                        continue;
                }
            }

            // === Per-dipole phase computation ===

            // Step 1: Back-project dipole to reflection facet along -seg.direction → R
            double dp_rn_refl = dip.x * nBx + dip.y * nBy + dip.z * nBz + nBd;
            if (fabs(dDotN) < 1e-10) continue;
            double alpha = -dp_rn_refl / dDotN;
            double Rx = dip.x + alpha * dDir.x;
            double Ry = dip.y + alpha * dDir.y;
            double Rz = dip.z + alpha * dDir.z;
            double dist_R_to_dip = fabs(alpha);

            // Step 2: Trace from R along -d_inc to find entry facet → P
            double minT_entry = 1e30;
            int entryFacetFound = -1;
            for (int f = 0; f < m_particle->nFacets; ++f)
            {
                if (f == seg.entryFacetId) continue;
                const Facet &fac = m_particle->facets[f];
                double nfx = (double)fac.in_normal.cx;
                double nfy = (double)fac.in_normal.cy;
                double nfz = (double)fac.in_normal.cz;
                double nfd = fac.in_normal.d_param;

                // Ray: R + t_f * (-d_inc), plane: n_f · x + d = 0
                // t_f = (n_f · R + d) / (n_f · d_inc)
                double dp_Rf = Rx * nfx + Ry * nfy + Rz * nfz + nfd;
                double dp_df = dIncX * nfx + dIncY * nfy + dIncZ * nfz;
                if (fabs(dp_df) < 1e-10) continue;

                double t_f = dp_Rf / dp_df;
                if (t_f < 1e-6 || t_f >= minT_entry) continue;

                minT_entry = t_f;
                entryFacetFound = f;
            }

            if (entryFacetFound < 0)
            {
                ++nEntryMissed;
                continue;
            }

            // P = R - minT_entry * d_inc
            double dist_PR = minT_entry;
            double Px = Rx - dist_PR * dIncX;
            double Py = Ry - dist_PR * dIncY;
            double Pz = Rz - dist_PR * dIncZ;

            // Step 3: Compute per-dipole OP
            // OP = FAR_ZONE + incDir·P + n * (dist_entry_to_refl + dist_refl_to_dip)
            double externalOP = FAR_ZONE_DISTANCE + iDirX * Px + iDirY * Py + iDirZ * Pz;
            double internalDist = dist_PR + dist_R_to_dip;
            double total_OP = externalOP + re_n * internalDist;

            complex phase = exp_im(m_k * total_OP);
            if (im_n > 1e-15)
            {
                phase = phase * exp(-m_k * im_n * internalDist);
            }

            // Diffraction weight: Kirchhoff integral for narrow beams,
            // product of half-planes for wide beams
            complex kirchWeight(1.0, 0.0);
            if (kirchhoffReady)
            {
                double Rdx = Rx - kb_ox, Rdy = Ry - kb_oy, Rdz = Rz - kb_oz;
                double obs_u = Rdx*e_perp.x + Rdy*e_perp.y + Rdz*e_perp.z;
                double obs_v = Rdx*e_par.x  + Rdy*e_par.y  + Rdz*e_par.z;

                double z_f = (dist_R_to_dip > 1e-10) ? dist_R_to_dip : 1e-10;
                double scale_f = sqrt(2.0 / (lambda_eff * z_f));

                // Max Fresnel-scaled distance from obs to polygon vertices
                double maxFresnel = 0;
                for (int vi = 0; vi < seg.polyNVertices; ++vi)
                {
                    double Xv = fabs((polyB_u[vi] - obs_u) * scale_f);
                    double Yv = fabs((polyB_v[vi] - obs_v) * scale_f);
                    if (Xv > maxFresnel) maxFresnel = Xv;
                    if (Yv > maxFresnel) maxFresnel = Yv;
                }

                if (maxFresnel > 2.5)
                {
                    // Large Fresnel number: product of half-plane weights (accurate)
                    kirchWeight = complex(1.0, 0.0);
                    for (int ei = 0; ei < seg.polyNVertices; ++ei)
                    {
                        int ej = (ei + 1) % seg.polyNVertices;
                        double eu = polyB_u[ej] - polyB_u[ei];
                        double ev = polyB_v[ej] - polyB_v[ei];
                        double elen = sqrt(eu*eu + ev*ev);
                        if (elen < 1e-12) continue;
                        double d_edge = ((obs_v - polyB_v[ei])*eu
                                       - (obs_u - polyB_u[ei])*ev) / elen;
                        kirchWeight = kirchWeight * FresnelEdgeWeight(d_edge, z_f, lambda_eff);
                    }
                }
                else
                {
                    // Small Fresnel number: exact Kirchhoff integral
                    kirchWeight = KirchhoffPolygonWeight(
                        polyB_u, polyB_v, seg.polyNVertices,
                        obs_u, obs_v, z_f, lambda_eff);
                }
            }
            if (norm(kirchWeight) < 1e-8)
                continue;

            // Additive accumulation: reflected wave adds to direct wave
            complex wPhase = kirchWeight * phase;

            // Y-polarization: pure s-input -> (0, 1) in (p, s) basis
            complex Ep_Y = seg.J.m22 * wPhase;
            complex Epar_Y = seg.J.m12 * wPhase;
            dip.Ex_Y += Ep_Y * e_perp.x + Epar_Y * e_par.x;
            dip.Ey_Y += Ep_Y * e_perp.y + Epar_Y * e_par.y;
            dip.Ez_Y += Ep_Y * e_perp.z + Epar_Y * e_par.z;

            // X-polarization: pure -p-input -> (-1, 0) in (p, s) basis
            complex Ep_X = -seg.J.m21 * wPhase;
            complex Epar_X = -seg.J.m11 * wPhase;
            dip.Ex_X += Ep_X * e_perp.x + Epar_X * e_par.x;
            dip.Ey_X += Ep_X * e_perp.y + Epar_X * e_par.y;
            dip.Ez_X += Ep_X * e_perp.z + Epar_X * e_par.z;

            if (norm(kirchWeight) < 0.99) ++nSmoothed;
            ++segUpdated;
        }
        nDipolesUpdated += segUpdated;
    }

    cout << "GO reflected beams (nActs=1): " << nProcessed << " segments processed, "
         << nSkippedNorm << " skipped (|J|>" << maxJonesNorm << "), "
         << nDipolesUpdated << " dipole updates, "
         << nSmoothed << " Kirchhoff-weighted, "
         << nEntryMissed << " entry-facet misses" << endl;
}

void ADDAFieldComputer::ApplyUniformBeam(const Point3f &refractedDir,
                                          const Point3f &incidentDir,
                                          double Ts, double Tp)
{
    // Apply a uniform refracted plane wave to ALL dipoles for non-convex particles.
    //
    // Phase: total_OP = FAR_ZONE + dot(incDir, r_ref) + n * dot(d_refr, r - r_ref)
    // where r_ref = (0, 0, m_zMax)
    // Ts, Tp: Fresnel transmission coefficients (computed from cap normal in caller)

    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    double dx = (double)refractedDir.cx;
    double dy = (double)refractedDir.cy;
    double dz = (double)refractedDir.cz;

    // Build s/p polarization basis from plane of incidence (incDir, d_refr)
    // s = normalize(cross(incDir, d_refr)) — perpendicular to plane of incidence
    double ix = (double)incidentDir.cx;
    double iy = (double)incidentDir.cy;
    double iz = (double)incidentDir.cz;

    double sx = iy * dz - iz * dy;
    double sy = iz * dx - ix * dz;
    double sz = ix * dy - iy * dx;
    double slen = sqrt(sx*sx + sy*sy + sz*sz);

    // Precompute E-field vectors for Y-pol and X-pol
    // These are constant across all dipoles (uniform plane wave)
    double EY_x, EY_y, EY_z;  // refracted E for Y-pol incident (0,1,0)
    double EX_x, EX_y, EX_z;  // refracted E for X-pol incident (-1,0,0)

    if (slen < 1e-6)
    {
        // Near-normal incidence: Ts ≈ Tp, no rotation needed
        double T = 0.5 * (Ts + Tp);
        EY_x = 0; EY_y = T; EY_z = 0;
        EX_x = -T; EX_y = 0; EX_z = 0;
    }
    else
    {
        sx /= slen; sy /= slen; sz /= slen;

        // p_inc = cross(s, incDir) — p-direction for incident beam
        double pix = sy * iz - sz * iy;
        double piy = sz * ix - sx * iz;
        double piz = sx * iy - sy * ix;
        double pilen = sqrt(pix*pix + piy*piy + piz*piz);
        pix /= pilen; piy /= pilen; piz /= pilen;

        // p_refr = cross(s, d_refr) — p-direction for refracted beam
        double prx = sy * dz - sz * dy;
        double pry = sz * dx - sx * dz;
        double prz = sx * dy - sy * dx;
        double prlen = sqrt(prx*prx + pry*pry + prz*prz);
        prx /= prlen; pry /= prlen; prz /= prlen;

        // Y-pol: E_inc = (0,1,0). Decompose into s and p_inc:
        double Eys = sy;       // dot((0,1,0), s)
        double Eyp = piy;     // dot((0,1,0), p_inc)
        EY_x = Ts * Eys * sx + Tp * Eyp * prx;
        EY_y = Ts * Eys * sy + Tp * Eyp * pry;
        EY_z = Ts * Eys * sz + Tp * Eyp * prz;

        // X-pol: E_inc = (-1,0,0). Decompose into s and p_inc:
        double Exs = -sx;     // dot((-1,0,0), s)
        double Exp = -pix;    // dot((-1,0,0), p_inc)
        EX_x = Ts * Exs * sx + Tp * Exp * prx;
        EX_y = Ts * Exs * sy + Tp * Exp * pry;
        EX_z = Ts * Exs * sz + Tp * Exp * prz;
    }

    static const double FAR_ZONE_DISTANCE = 10000.0;
    double externalOP = FAR_ZONE_DISTANCE + iz * m_zMax;
    double dRefDotRef = dz * m_zMax;

    for (size_t di = 0; di < m_dipoles.size(); ++di)
    {
        DipoleField &dip = m_dipoles[di];
        double t = dx * dip.x + dy * dip.y + dz * dip.z - dRefDotRef;

        double total_OP = externalOP + re_n * t;
        complex phase = exp_im(m_k * total_OP);

        if (im_n > 1e-15)
        {
            double tAbs = (t > 0) ? t : 0;
            phase = phase * exp(-m_k * im_n * tAbs);
        }

        dip.nContributions = 1;

        dip.Ex_Y = EY_x * phase;
        dip.Ey_Y = EY_y * phase;
        dip.Ez_Y = EY_z * phase;

        dip.Ex_X = EX_x * phase;
        dip.Ey_X = EX_y * phase;
        dip.Ez_X = EX_z * phase;
    }

    cout << "Applied uniform beam to all " << m_dipoles.size() << " dipoles"
         << " (d_refr=" << dx << "," << dy << "," << dz
         << ", Ts=" << Ts << ", Tp=" << Tp << ")" << endl;
    cout << "  EY=(" << EY_x << "," << EY_y << "," << EY_z
         << ") EX=(" << EX_x << "," << EX_y << "," << EX_z << ")" << endl;
}

void ADDAFieldComputer::ResetAllFields()
{
    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        DipoleField &dip = m_dipoles[i];
        dip.Ex_Y = dip.Ey_Y = dip.Ez_Y = complex(0, 0);
        dip.Ex_X = dip.Ey_X = dip.Ez_X = complex(0, 0);
        dip.nContributions = 0;
    }
}

void ADDAFieldComputer::DiagnoseGOvsPW(const Point3f &incidentDir)
{
    // 1. Print contribution statistics
    int n0 = 0, n1 = 0, n2 = 0, n3plus = 0;
    for (size_t i = 0; i < m_dipoles.size(); ++i)
    {
        int nc = m_dipoles[i].nContributions;
        if (nc == 0) ++n0;
        else if (nc == 1) ++n1;
        else if (nc == 2) ++n2;
        else ++n3plus;
    }
    cout << "\n=== GO Field Diagnostics ===" << endl;
    cout << "Contributions: 0=" << n0 << " 1=" << n1
         << " 2=" << n2 << " 3+=" << n3plus << endl;

    // 2. For GO-covered dipoles, compute per-facet PW independently and compare
    double re_n = real(m_ri);
    double ix = (double)incidentDir.cx;
    double iy = (double)incidentDir.cy;
    double iz = (double)incidentDir.cz;

    static const double FAR_ZONE_DISTANCE = 10000.0;

    // Build per-facet refraction data (same as FillUncoveredPerFacet)
    struct FacetInfo {
        int id;
        double cosI;
        Point3f refDir;
        double Ts, Tp;
        double EY_x, EY_y, EY_z;
        double EX_x, EX_y, EX_z;
        double phaseRef;
    };
    vector<FacetInfo> facets;

    for (int f = 0; f < m_particle->nFacets; ++f)
    {
        const Facet &fac = m_particle->facets[f];
        double nx = (double)fac.in_normal.cx;
        double ny = (double)fac.in_normal.cy;
        double nz = (double)fac.in_normal.cz;
        double cosI = ix*nx + iy*ny + iz*nz;
        if (cosI < 1e-6) continue;

        double sinIsq = 1.0 - cosI*cosI;
        double cosT = sqrt(fmax(0.0, 1.0 - sinIsq/(re_n*re_n)));
        double eta = 1.0 / re_n;
        double C = cosT - eta * cosI;
        double rdx = eta*ix + C*nx;
        double rdy = eta*iy + C*ny;
        double rdz = eta*iz + C*nz;
        double rlen = sqrt(rdx*rdx + rdy*rdy + rdz*rdz);
        rdx /= rlen; rdy /= rlen; rdz /= rlen;

        double Ts = 2.0*cosI / (cosI + re_n*cosT);
        double Tp = 2.0*cosI / (re_n*cosI + cosT);

        double sx = iy*rdz - iz*rdy;
        double sy = iz*rdx - ix*rdz;
        double sz = ix*rdy - iy*rdx;
        double slen = sqrt(sx*sx + sy*sy + sz*sz);

        double EY_x, EY_y, EY_z;
        double EX_x, EX_y, EX_z;

        if (slen < 1e-6)
        {
            double T = 0.5*(Ts+Tp);
            EY_x=0; EY_y=T; EY_z=0;
            EX_x=-T; EX_y=0; EX_z=0;
        }
        else
        {
            sx/=slen; sy/=slen; sz/=slen;
            double pix_=sy*iz-sz*iy, piy_=sz*ix-sx*iz, piz_=sx*iy-sy*ix;
            double pilen=sqrt(pix_*pix_+piy_*piy_+piz_*piz_);
            pix_/=pilen; piy_/=pilen; piz_/=pilen;
            double prx=sy*rdz-sz*rdy, pry=sz*rdx-sx*rdz, prz=sx*rdy-sy*rdx;
            double prlen=sqrt(prx*prx+pry*pry+prz*prz);
            prx/=prlen; pry/=prlen; prz/=prlen;

            double Eys=sy, Eyp=piy_;
            EY_x=Ts*Eys*sx+Tp*Eyp*prx;
            EY_y=Ts*Eys*sy+Tp*Eyp*pry;
            EY_z=Ts*Eys*sz+Tp*Eyp*prz;
            double Exs=-sx, Exp_=-pix_;
            EX_x=Ts*Exs*sx+Tp*Exp_*prx;
            EX_y=Ts*Exs*sy+Tp*Exp_*pry;
            EX_z=Ts*Exs*sz+Tp*Exp_*prz;
        }

        FacetInfo fi;
        fi.id=f; fi.cosI=cosI;
        fi.refDir=Point3f((float)rdx,(float)rdy,(float)rdz);
        fi.Ts=Ts; fi.Tp=Tp;
        fi.EY_x=EY_x; fi.EY_y=EY_y; fi.EY_z=EY_z;
        fi.EX_x=EX_x; fi.EX_y=EX_y; fi.EX_z=EX_z;
        double cx_=(double)fac.center.cx, cy_=(double)fac.center.cy, cz_=(double)fac.center.cz;
        fi.phaseRef=FAR_ZONE_DISTANCE+ix*cx_+iy*cy_+iz*cz_;
        facets.push_back(fi);
    }

    // For the first 10 GO-covered dipoles, compare GO vs PW
    int nPrinted = 0;
    double sumErr2_Y = 0, sumNorm2_Y = 0;
    double sumErr2_X = 0, sumNorm2_X = 0;
    int nCompared = 0;

    for (size_t i = 0; i < m_dipoles.size() && nCompared < 5000; ++i)
    {
        DipoleField &dip = m_dipoles[i];
        if (dip.nContributions == 0) continue;

        // Find the matching facet for this dipole (back-project)
        bool foundFacet = false;
        for (size_t fi = 0; fi < facets.size(); ++fi)
        {
            const FacetInfo &fr = facets[fi];
            const Facet &fac = m_particle->facets[fr.id];
            double dx=(double)fr.refDir.cx, dy=(double)fr.refDir.cy, dz=(double)fr.refDir.cz;

            double dp_nr = dip.x*fac.in_normal.cx + dip.y*fac.in_normal.cy
                         + dip.z*fac.in_normal.cz + fac.in_normal.d_param;
            double dp_nd = dx*fac.in_normal.cx + dy*fac.in_normal.cy + dz*fac.in_normal.cz;
            if (fabs(dp_nd) < 1e-10) continue;
            double t_back = -dp_nr / dp_nd;
            if (t_back > 1e-6) continue;

            Point3f proj((float)(dip.x+dx*t_back),(float)(dip.y+dy*t_back),(float)(dip.z+dz*t_back));
            Point3f polyNormal(0,0,0);
            if (fac.nVertices >= 3)
            {
                Point3f e1=fac.arr[1]-fac.arr[0], e2=fac.arr[2]-fac.arr[0];
                polyNormal=CrossProduct(e1,e2);
                float pnLen=(float)Length(polyNormal);
                if (pnLen > 1e-10f) polyNormal=polyNormal/(double)pnLen;
            }
            if (!IsInsidePolygon(proj, fac.arr, fac.nVertices, polyNormal)) continue;

            // Compute PW field for this dipole
            double cx_=(double)fac.center.cx, cy_=(double)fac.center.cy, cz_=(double)fac.center.cz;
            double t_int = dx*(dip.x-cx_) + dy*(dip.y-cy_) + dz*(dip.z-cz_);
            double total_OP = fr.phaseRef + re_n * t_int;
            complex phase = exp_im(m_k * total_OP);

            complex pw_Ex_Y = fr.EY_x * phase;
            complex pw_Ey_Y = fr.EY_y * phase;
            complex pw_Ez_Y = fr.EY_z * phase;
            complex pw_Ex_X = fr.EX_x * phase;
            complex pw_Ey_X = fr.EX_y * phase;
            complex pw_Ez_X = fr.EX_z * phase;

            // Accumulate error
            sumErr2_Y += norm(dip.Ex_Y - pw_Ex_Y) + norm(dip.Ey_Y - pw_Ey_Y) + norm(dip.Ez_Y - pw_Ez_Y);
            sumNorm2_Y += norm(pw_Ex_Y) + norm(pw_Ey_Y) + norm(pw_Ez_Y);
            sumErr2_X += norm(dip.Ex_X - pw_Ex_X) + norm(dip.Ey_X - pw_Ey_X) + norm(dip.Ez_X - pw_Ez_X);
            sumNorm2_X += norm(pw_Ex_X) + norm(pw_Ey_X) + norm(pw_Ez_X);
            ++nCompared;

            // Print first few
            if (nPrinted < 5)
            {
                cout << "  Dip#" << i << " nc=" << dip.nContributions
                     << " r=(" << dip.x << "," << dip.y << "," << dip.z
                     << ") facet=" << fr.id << endl;
                cout << "    GO  Y: (" << real(dip.Ex_Y) << "," << real(dip.Ey_Y) << "," << real(dip.Ez_Y)
                     << ") + i(" << imag(dip.Ex_Y) << "," << imag(dip.Ey_Y) << "," << imag(dip.Ez_Y) << ")" << endl;
                cout << "    PW  Y: (" << real(pw_Ex_Y) << "," << real(pw_Ey_Y) << "," << real(pw_Ez_Y)
                     << ") + i(" << imag(pw_Ex_Y) << "," << imag(pw_Ey_Y) << "," << imag(pw_Ez_Y) << ")" << endl;
                cout << "    GO  X: (" << real(dip.Ex_X) << "," << real(dip.Ey_X) << "," << real(dip.Ez_X)
                     << ") + i(" << imag(dip.Ex_X) << "," << imag(dip.Ey_X) << "," << imag(dip.Ez_X) << ")" << endl;
                cout << "    PW  X: (" << real(pw_Ex_X) << "," << real(pw_Ey_X) << "," << real(pw_Ez_X)
                     << ") + i(" << imag(pw_Ex_X) << "," << imag(pw_Ey_X) << "," << imag(pw_Ez_X) << ")" << endl;
                ++nPrinted;
            }
            foundFacet = true;
            break;
        }
    }

    if (nCompared > 0)
    {
        double relErr_Y = sqrt(sumErr2_Y / sumNorm2_Y);
        double relErr_X = sqrt(sumErr2_X / sumNorm2_X);
        cout << "GO vs PW comparison (" << nCompared << " dipoles):" << endl;
        cout << "  Y-pol relative error: " << relErr_Y << endl;
        cout << "  X-pol relative error: " << relErr_X << endl;
    }
    cout << "=== End Diagnostics ===\n" << endl;
}
