#pragma once

#include "geometry_lib.h"
#include "Handler.h"

class Tracer
{
public:
    Tracer(Particle *particle, int nActs, const std::string &resultFileName);
    ~Tracer();

    void SetHandler(Handler *handler);

    Light m_incidentLight;
    std::string m_summary;
    Scattering *m_scattering;

protected:
    Handler *m_handler;
    Particle *m_particle;

    std::string m_resultDirName;
    Symmetry m_symmetry;

private:
    void SetIncidentLight(Particle *particle);
};
