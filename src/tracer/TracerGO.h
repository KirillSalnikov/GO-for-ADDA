#pragma once

#include "Tracer.h"

class TracerGO : public Tracer
{
public:
	TracerGO(Particle *particle, int reflNum, const std::string &resultFileName);

	void TraceFixed(const double &beta, const double &gamma);
};
