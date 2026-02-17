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

    int nActs;            // number of internal reflections (beam generation)
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

    int nContributions;  // number of GO beam segments covering this dipole
    int assignedFacet;   // entry facet assigned by FillUncoveredPerFacet (-1 = none)

    DipoleField()
        : x(0), y(0), z(0), ix(0), iy(0), iz(0),
          Ex_Y(0,0), Ey_Y(0,0), Ez_Y(0,0),
          Ex_X(0,0), Ey_X(0,0), Ez_X(0,0),
          nContributions(0), assignedFacet(-1)
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

    /// Fill uncovered dipoles (nContributions==0) with refracted plane wave
    void FillUncoveredWithPlaneWave(const Point3f &refractedDir,
                                     const Point3f &incidentDir,
                                     double Ts, double Tp);

    /// Fill uncovered dipoles with per-facet refracted plane waves.
    /// Each illuminated facet gets its own Snell refraction + Fresnel coefficients.
    void FillUncoveredPerFacet(const Point3f &incidentDir);

    /// Apply a uniform refracted plane wave to ALL dipoles (for non-convex).
    /// Uses the given refracted direction with z_max-based phase and
    /// analytical Fresnel transmission with full 3D polarization.
    void ApplyUniformBeam(const Point3f &refractedDir,
                          const Point3f &incidentDir,
                          double Ts, double Tp);

    int GetNumDipoles() const { return (int)m_dipoles.size(); }
    int GetNumSegments() const { return m_segCount; }
    double GetGridSpacing() const { return m_gridSpacing; }

    void IncrementSegCount() { ++m_segCount; }

    /// Add per-facet first reflection (Fabry-Perot correction).
    /// For each dipole assigned to an entry facet, traces refracted beam to exit facet,
    /// computes Fresnel reflection, and adds reflected wave. No overlap.
    void AddPerFacetReflection(const Point3f &incidentDir);

    /// Accumulate E-field from GO-traced nActs=1 reflected beam segments.
    /// Uses per-dipole phase: back-projects each dipole to reflection facet,
    /// then to entry facet, computing OP = FAR_ZONE + incDir·P + n*(d1+d2).
    void AccumulateReflectedBeams(const std::vector<InternalBeamSegment> &segments,
                                   const Point3f &incidentDir,
                                   double maxJonesNorm = 0.5,
                                   bool useDiffraction = false,
                                   double minFresnelNum = 1.0);

    /// Print contribution stats and compare GO field vs per-facet PW for covered dipoles
    void DiagnoseGOvsPW(const Point3f &incidentDir);

    /// Reset all dipole fields (for diagnostic re-computation)
    void ResetAllFields();

private:
    bool IsInsideParticle(const Point3f &p) const;
    bool IsInsidePolygon(const Point3f &proj, const Point3f *polyArr,
                         int nVert, const Point3f &normal) const;
    /// Minimum distance from point to nearest polygon edge (for boundary smoothing)
    double MinEdgeDistance(const Point3f &proj, const Point3f *polyArr, int nVert) const;

    Particle *m_particle;
    double m_wavelength;
    complex m_ri;
    double m_gridSpacing;   // dipole spacing = wavelength / dpl
    double m_k;             // wave number = 2*PI / wavelength

    std::vector<DipoleField> m_dipoles;
    int m_segCount;
    double m_zMax;         // Top of particle (for consistent phase reference)

    // Grid bounding box dimensions
    int m_boxX, m_boxY, m_boxZ;
};
