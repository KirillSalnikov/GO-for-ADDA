#include "ADDAField.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

ADDAFieldComputer::ADDAFieldComputer(Particle *particle, double wavelength,
                                     const complex &ri, int dpl)
    : m_particle(particle),
      m_wavelength(wavelength),
      m_ri(ri),
      m_segCount(0),
      m_boxX(0), m_boxY(0), m_boxZ(0)
{
    m_gridSpacing = wavelength / dpl;
    m_k = 2.0 * M_PI / wavelength;
}

bool ADDAFieldComputer::IsInsideParticle(const Point3f &p) const
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

void ADDAFieldComputer::AccumulateBeamContribution(const InternalBeamSegment &seg)
{
    // Beam local polarization basis
    Point3d e_perp(seg.polarizationBasis.cx,
                   seg.polarizationBasis.cy,
                   seg.polarizationBasis.cz);
    Point3d d(seg.direction.cx, seg.direction.cy, seg.direction.cz);
    Point3d e_par = CrossProductD(e_perp, d);
    double len = LengthD(e_par);
    if (len < 1e-12) return;
    e_par = e_par / len;

    double re_n = real(m_ri);
    double im_n = imag(m_ri);

    // Compute polygon normal from vertices (determines winding direction)
    Point3f polyNormal(0, 0, 0);
    if (seg.polyNVertices >= 3)
    {
        Point3f e1 = seg.polyArr[1] - seg.polyArr[0];
        Point3f e2 = seg.polyArr[2] - seg.polyArr[0];
        polyNormal = CrossProduct(e1, e2);
        float pnLen = (float)Length(polyNormal);
        if (pnLen > 1e-10f)
        {
            polyNormal = polyNormal / (double)pnLen;
        }
    }

    // Exit facet plane for projection
    Point3f exitNorm = seg.exitNormal;

    for (size_t di = 0; di < m_dipoles.size(); ++di)
    {
        DipoleField &dip = m_dipoles[di];
        Point3f r((float)dip.x, (float)dip.y, (float)dip.z);

        // Distance from beam entry plane to dipole along beam direction
        // beam.front = DotProduct(-direction, beamCenter)
        // t = DotProduct(direction, r) + front = DotProduct(direction, r - beamCenter)
        double t = (double)seg.direction.cx * r.cx
                 + (double)seg.direction.cy * r.cy
                 + (double)seg.direction.cz * r.cz
                 + seg.front;

        // Check within segment bounds
        if (t < -1e-6 || t > seg.segmentLength + 1e-6)
            continue;

        // Project dipole to exit facet plane along beam direction
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
        // Use polygon's own normal (from vertex winding) for correct containment test
        if (!IsInsidePolygon(r_proj, seg.polyArr, seg.polyNVertices, polyNormal))
            continue;

        // Compute phase
        double delta_OP = t * re_n;
        double total_OP = seg.opticalPath + delta_OP;
        complex phase = exp_im(m_k * total_OP);

        // Absorption
        if (im_n > 1e-15)
        {
            phase = phase * exp(-m_k * im_n * t);
        }

        // Y-polarization: E_in = (1, 0) in (e_perp, e_par) basis
        complex Ep_Y = seg.J.m11 * phase;   // e_perp component
        complex Epar_Y = seg.J.m21 * phase;  // e_par component
        dip.Ex_Y += Ep_Y * e_perp.x + Epar_Y * e_par.x;
        dip.Ey_Y += Ep_Y * e_perp.y + Epar_Y * e_par.y;
        dip.Ez_Y += Ep_Y * e_perp.z + Epar_Y * e_par.z;

        // X-polarization: E_in = (0, +1) in (e_perp, e_par) basis
        // ADDA's incPolX = (-1,0,0) when prop=(0,0,-1), and e_par = -X, so E·e_par = +1
        complex Ep_X = seg.J.m12 * phase;
        complex Epar_X = seg.J.m22 * phase;
        dip.Ex_X += Ep_X * e_perp.x + Epar_X * e_par.x;
        dip.Ey_X += Ep_X * e_perp.y + Epar_X * e_par.y;
        dip.Ez_X += Ep_X * e_perp.z + Epar_X * e_par.z;
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
