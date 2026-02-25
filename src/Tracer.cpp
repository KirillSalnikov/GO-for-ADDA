#include "Tracer.h"

#include <string>
#include <iostream>

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

using namespace std;

Tracer::Tracer(Particle *particle, int nActs, const string &resultFileName)
    : m_resultDirName(resultFileName)
{
    SetIncidentLight(particle);

    if (particle->IsConcave())
    {
        m_scattering = new ScatteringNonConvex(particle, &m_incidentLight, true, nActs);
    }
    else
    {
        m_scattering = new ScatteringConvex(particle, &m_incidentLight, true, nActs);
    }

    m_particle = m_scattering->m_particle;
    m_symmetry = m_particle->GetSymmetry();
}

Tracer::~Tracer()
{
}

void Tracer::SetIncidentLight(Particle *particle)
{
    m_incidentLight.direction = Point3f(0, 0, -1);
    m_incidentLight.polarizationBasis = Point3f(0, 1, 0);

    Point3f point = m_incidentLight.direction * particle->LongRadius();
    m_incidentLight.direction.d_param = DotProduct(point, m_incidentLight.direction);
}

void Tracer::SetHandler(Handler *handler)
{
    m_handler = handler;
    m_handler->SetScattering(m_scattering);
}
