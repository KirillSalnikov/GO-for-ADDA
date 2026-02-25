#include "TracerGO.h"
#include "HandlerGO.h"
#include <iostream>

using namespace std;

TracerGO::TracerGO(Particle *particle, int reflNum, const std::string &resultFileName)
	: Tracer(particle, reflNum, resultFileName)
{
}

void TracerGO::TraceFixed(const double &beta, const double &gamma)
{
	double b = DegToRad(beta);
	double g = DegToRad(gamma);

	vector<Beam> outBeams;
	m_scattering->ScatterLight(b, g, outBeams);
    m_handler->HandleBeams(outBeams, 0);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

    m_handler->WriteMatricesToFile(m_resultDirName, 1000);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

