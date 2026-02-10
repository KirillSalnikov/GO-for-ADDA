#pragma once

#include "geometry_lib.h"
#include "compl.hpp"
#include "Beam.h"
#include "Particle.h"
#include <vector>
#include <string>

/// Internal beam segment between two facet interactions
struct InternalBeamSegment
{
    // Beam state at entry (before Fresnel modification at exit facet)
    Matrix2x2c J;
    Point3f direction;
    Point3f polarizationBasis;
    double opticalPath;
    double front;

    // Intersection polygon on exit facet (beam cross-section)
    Point3f polyArr[MAX_VERTEX_NUM];
    int polyNVertices;

    // Facet info
    int entryFacetId;
    int exitFacetId;
    double segmentLength;

    // Exit facet plane for projection
    Point3f exitNormal;   // in_normal of exit facet
    double exitD;         // d_param of exit facet in_normal
};

/// Electric field accumulated at a single dipole
struct DipoleField
{
    double x, y, z;         // physical coordinates
    int ix, iy, iz;         // integer grid indices

    // Accumulated complex E-field for Y-polarization
    complex Ex_Y, Ey_Y, Ez_Y;
    // Accumulated complex E-field for X-polarization
    complex Ex_X, Ey_X, Ez_X;

    DipoleField()
        : x(0), y(0), z(0), ix(0), iy(0), iz(0),
          Ex_Y(0,0), Ey_Y(0,0), Ez_Y(0,0),
          Ex_X(0,0), Ey_X(0,0), Ez_X(0,0)
    {}
};

class ADDAFieldComputer
{
public:
    ADDAFieldComputer(Particle *particle, double wavelength,
                      const complex &ri, int dpl);

    /// Build cubic dipole grid inside convex particle
    void BuildDipoleGrid();

    /// Accumulate E-field from one internal beam segment to all covered dipoles
    void AccumulateBeamContribution(const InternalBeamSegment &seg);

    /// Write DDSCAT6-format geometry file for ADDA (-shape read)
    void WriteGeometryFile(const std::string &filename) const;

    /// Write init field file for Y-polarization (-init_field read)
    void WriteFieldFileY(const std::string &filename) const;

    /// Write init field file for X-polarization
    void WriteFieldFileX(const std::string &filename) const;

    int GetNumDipoles() const { return (int)m_dipoles.size(); }
    int GetNumSegments() const { return m_segCount; }
    double GetGridSpacing() const { return m_gridSpacing; }

    void IncrementSegCount() { ++m_segCount; }

private:
    bool IsInsideParticle(const Point3f &p) const;
    bool IsInsidePolygon(const Point3f &proj, const Point3f *polyArr,
                         int nVert, const Point3f &normal) const;

    Particle *m_particle;
    double m_wavelength;
    complex m_ri;
    double m_gridSpacing;   // dipole spacing = wavelength / dpl
    double m_k;             // wave number = 2*PI / wavelength

    std::vector<DipoleField> m_dipoles;
    int m_segCount;

    // Grid bounding box dimensions
    int m_boxX, m_boxY, m_boxZ;
};
