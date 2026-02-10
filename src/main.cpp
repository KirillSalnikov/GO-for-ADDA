#include <iostream>
#include <assert.h>
#include <float.h>
#include <string>

#include "HexagonalAggregate.h"
#include "ConcaveHexagonal.h"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "global.h"
#include "TracerGO.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"
#include "Droxtal.h"
#include "ADDAField.h"
#include "ScatteringConvex.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

enum class ParticleType : int
{
    Hexagonal = 1,
    Bullet = 2,
    BulletRosette = 3,
    Droxtal = 4,
    ConcaveHexagonal = 10,
    TiltedHexagonal = 11,
    HexagonalAggregate = 12,
    CertainAggregate = 999
};

Tracks trackGroups;

void SetArgRules(ArgPP &parser)
{
    int zero = 0;
    parser.AddRule("p", '+'); // particle (type, size, ...)
    parser.AddRule("ri", 2); // refractive index (Re and Im parts)
    parser.AddRule("n", 1); // number of internal reflection
    parser.AddRule("pf", zero, true); // particle (filename)
    parser.AddRule("rs", 1, true, "pf"); // resize particle (new size)
    parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
    parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
    parser.AddRule("w", 1, true); // wavelength
    parser.AddRule("grid", '+', true); /* backscattering grid:
 (radius, Nphi, Ntheta) when 3 parameters
 (theta1, theta2, Nphi, Ntheta) when 4 parameters*/
    parser.AddRule("tr", 1, true); // file with trajectories
    parser.AddRule("all", 0, true); // calculate all trajectories
    parser.AddRule("abs", zero, true, "w"); // accounting of absorbtion
    parser.AddRule("close", 0, true); // closing of program after calculation
    parser.AddRule("o", 1, true); // output folder name
    parser.AddRule("gr", zero, true);
    parser.AddRule("incoh", zero, true);
    parser.AddRule("forced_nonconvex", zero, true);
    parser.AddRule("forced_convex", zero, true);
    parser.AddRule("r", 1, true); // restriction ratio for small beams when intersection (100 by default)
    parser.AddRule("log", 1, true); // time of writing progress (in seconds)
    parser.AddRule("adda", 0, true); // ADDA mode: compute internal field and output for ADDA
    parser.AddRule("dpl", 1, true); // dipoles per wavelength for ADDA grid (default 10)
}

ScatteringRange SetConus(ArgPP &parser)
{
    if (parser.GetArgNumber("grid") == 3)
    {
        double radius = parser.GetDoubleValue("grid", 0);
        int nAz = parser.GetDoubleValue("grid", 1);
        int nZen = parser.GetDoubleValue("grid", 2);
        return ScatteringRange(M_PI - DegToRad(radius), M_PI, nAz, nZen);
    }
    else if (parser.GetArgNumber("grid") == 4)
    {
        double zenStart = parser.GetDoubleValue("grid", 0);
        double zenEnd = parser.GetDoubleValue("grid", 1);

        if (zenStart < zenEnd)
        {
            int nAz = parser.GetDoubleValue("grid", 2);
            int nZen = parser.GetDoubleValue("grid", 3);
            return ScatteringRange(DegToRad(zenStart), DegToRad(zenEnd), nAz, nZen);
        }
        else
        {
            std::cerr << "ERROR!!!!!!! In \"grid\" arg 1 must be less than arg 2";
            exit(1);
        }
    }
    else
    {
        std::cerr << "ERROR!!!!!!! Wrong \"grid\" argument number";
        exit(1);
    }
}

AngleRange GetRange(const ArgPP &parser, const std::string &key,
                    Particle *particle)
{
    int number;
    double min, max;

    if (key == "b")
    {
        number = parser.GetIntValue("random", 0);
        min = 0;
        max = particle->GetSymmetry().beta;
    }
    else if (key == "g")
    {
        number = parser.GetIntValue("random", 1);
        min = 0;
        max = particle->GetSymmetry().gamma;
    }
    else
    {
        cerr << "Error! " << __FUNCTION__;
        throw std::exception();
    }

    return AngleRange(min, max, number);
}

int main(int argc, const char* argv[])
{
    RenameConsole("MBS");

    std::string additionalSummary;

    if (argc <= 1) // no arguments
    {
        cout << "No arguments. Press any key to exit..." << endl;
        getchar();
        return 1;
    }

    ArgPP args;
    SetArgRules(args);
    args.Parse(argc, argv);

    bool isAbs = args.IsCatched("abs");

    double re = args.GetDoubleValue("ri", 0);
    double im = args.GetDoubleValue("ri", 1);
    complex refrIndex = complex(re, im);

    if (!isAbs)
    {
        refrIndex = complex(re, 0);
    }

    Particle *particle = nullptr;


    additionalSummary += "Particle: ";

    if (args.IsCatched("pf"))
    {
        std::string filename = args.GetStringValue("p");
        particle = new Particle();
        particle->SetFromFile(filename);
        particle->SetRefractiveIndex(complex(refrIndex));

        double origDMax = particle->MaximalDimention();
        additionalSummary += "from file: " + filename + "\n";
        additionalSummary += "\tOriginal Dmax: " + std::to_string(origDMax);

        if (args.IsCatched("rs"))
        {
            double newSize = args.GetDoubleValue("rs", 0);
            particle->Resize(newSize);
        }

        double newDMax = particle->MaximalDimention();
        additionalSummary += ", new Dmax: " + std::to_string(newDMax)
                + ", resize factor: " + std::to_string(newDMax/origDMax) + '\n';
    }
    else
    {
        ParticleType type = (ParticleType)args.GetIntValue("p", 0);
        double height;
        double diameter;

        double sup;
        int num;

        switch (type)
        {
        case ParticleType::Hexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            particle = new Hexagonal(refrIndex, diameter, height);
            break;
        case ParticleType::Bullet:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (diameter*sqrt(3)*tan(DegToRad(62)))/4;
            particle = new Bullet(refrIndex, diameter, height, sup);
            break;
        case ParticleType::BulletRosette:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = (args.GetArgNumber("p") == 4)
                    ? args.GetDoubleValue("p", 3)
                    : diameter*sqrt(3)*tan(DegToRad(62))/4;
            particle = new BulletRosette(refrIndex, diameter, height, sup);
            break;
        case ParticleType::Droxtal:
            sup = args.GetDoubleValue("p", 3);
            particle = new Droxtal(refrIndex, DegToRad(32.35), DegToRad(71.81), sup);
            break;
        case ParticleType::ConcaveHexagonal:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            sup = args.GetDoubleValue("p", 3);
            particle = new ConcaveHexagonal(refrIndex, diameter, height, sup);
            break;
        case ParticleType::HexagonalAggregate:
            height = args.GetDoubleValue("p", 1);
            diameter = args.GetDoubleValue("p", 2);
            num = args.GetIntValue("p", 3);
            particle = new HexagonalAggregate(refrIndex, diameter, height, num);
            break;
        case ParticleType::CertainAggregate:
            sup = args.GetDoubleValue("p", 1);
            particle = new CertainAggregate(refrIndex, sup);
            break;
        default:
            assert(false && "ERROR! Incorrect type of particle.");
            break;
        }
    }

    additionalSummary += "\tRefractive index: " + to_string(re);

    if (fabs(im) > FLT_EPSILON)
    {
        additionalSummary +=  " + i\n";
    }

    additionalSummary += "\n";

    particle->Output("particle_for_check.dat");
    additionalSummary += "\tArea:" + to_string(particle->Area()) + "\n\n";

    int reflNum = args.GetDoubleValue("n");
    additionalSummary += "Number of secondary reflections: " + to_string(reflNum) + "\n";

    string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
                                           : "M";
    size_t pos = 0;

    while ((pos = dirName.find('%', pos)) != string::npos)
    {
        size_t start = pos;
        ++pos;
        string key;

        while (dirName[pos] != '_' && pos < dirName.size())
        {
            key += dirName[pos++];
        }

        if (key.size())
        {
            string val = args.GetStringValue(key);
            dirName.replace(start, pos-start, val);
            pos = start + val.size();
        }
    }

    bool isOutputGroups = args.IsCatched("gr");
    double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
    additionalSummary += "Wavelength (um): " + to_string(wave) + "\n";

    if (args.IsCatched("forced_nonconvex"))
    {
        particle->isConcave = true;
    }

    if (args.IsCatched("forced_convex"))
    {
        particle->isConcave = false;
    }

    if (args.IsCatched("tr"))
    {
        string trackFileName = args.GetStringValue("tr");
        trackGroups.ImportTracks(particle->nFacets, trackFileName);
        trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
        trackGroups.shouldOutputGroups = args.IsCatched("gr");
    }

    int nTheta = args.GetDoubleValue("grid", 2);

    additionalSummary += "Method: Geometrical optics";

    TracerGO tracer(particle, reflNum, dirName);
    tracer.m_scattering->m_wave = wave;
    if (args.IsCatched("r"))
    {
        tracer.m_scattering->restriction = args.GetDoubleValue("r", 0);
    }
    tracer.m_summary = additionalSummary;
    tracer.SetIsOutputGroups(isOutputGroups);

    HandlerGO *handler;

    if (args.IsCatched("tr"))
    {
        handler = new HandlerTracksGO(particle, &tracer.m_incidentLight, nTheta, wave);
        handler->SetTracks(&trackGroups);
    }
    else
    {
        handler = new HandlerTotalGO(particle, &tracer.m_incidentLight, nTheta, wave);
    }

    tracer.m_summary = additionalSummary;
    ScatteringRange grid = SetConus(args);
    handler->isCoh = !args.IsCatched("incoh");
    handler->SetScatteringSphere(grid);
    handler->SetAbsorptionAccounting(isAbs);
    tracer.SetHandler(handler);

    bool addaMode = args.IsCatched("adda");

    if (addaMode)
    {
        if (!args.IsCatched("fixed"))
        {
            cerr << "ERROR: -adda requires -fixed orientation" << endl;
            exit(1);
        }
        if (!args.IsCatched("w") || wave < 1e-15)
        {
            cerr << "ERROR: -adda requires -w (wavelength)" << endl;
            exit(1);
        }
    }

    if (args.IsCatched("fixed"))
    {
        additionalSummary += ", fixed orientation, ";
        double beta  = args.GetDoubleValue("fixed", 0);
        double gamma = args.GetDoubleValue("fixed", 1);
        additionalSummary += "zenith " + to_string(beta) + "\xC2\xB0, azimuth " + to_string(gamma) + "\xC2\xB0" + "\n\n";
        cout << additionalSummary;
        tracer.m_summary = additionalSummary;

        if (addaMode)
        {
            int dpl = args.IsCatched("dpl") ? args.GetIntValue("dpl") : 10;

            // Rotate particle to the specified orientation
            double b = DegToRad(beta);
            double g = DegToRad(gamma);
            particle->Rotate(b, g, 0);

            // Build dipole grid on the rotated particle
            ADDAFieldComputer addaField(particle, wave, refrIndex, dpl);
            addaField.BuildDipoleGrid();

            // Set up segment collection
            std::vector<InternalBeamSegment> segments;
            ScatteringConvex *sc = dynamic_cast<ScatteringConvex*>(tracer.m_scattering);
            if (sc)
            {
                sc->SetInternalSegmentStorage(&segments);
            }
            else
            {
                cerr << "ERROR: -adda only works with convex particles" << endl;
                exit(1);
            }

            // Run GO tracing (particle already rotated)
            vector<Beam> outBeams;
            tracer.m_scattering->ScatterLight(b, g, outBeams);
            handler->HandleBeams(outBeams, 0);
            outBeams.clear();
            handler->WriteMatricesToFile(dirName, 1000);

            cout << "Captured " << segments.size() << " internal beam segments" << endl;

            // Accumulate field contributions from all beam segments
            for (size_t si = 0; si < segments.size(); ++si)
            {
                addaField.AccumulateBeamContribution(segments[si]);
            }

            // Output ADDA files
            addaField.WriteGeometryFile(dirName + "_shape.dat");
            addaField.WriteFieldFileY(dirName + "_fieldY.dat");
            addaField.WriteFieldFileX(dirName + "_fieldX.dat");

            // Clean up
            sc->SetInternalSegmentStorage(nullptr);

            cout << "\nADDA usage:" << endl;
            cout << "  adda -shape read " << dirName << "_shape.dat"
                 << " -dpl " << dpl
                 << " -m " << re << " " << im
                 << " -prop 0 0 -1"
                 << " -init_field read " << dirName << "_fieldY.dat"
                 << " " << dirName << "_fieldX.dat"
                 << endl;
        }
        else
        {
            tracer.TraceFixed(beta, gamma);
        }
    }
    else if (args.IsCatched("random"))
    {
        if (addaMode)
        {
            cerr << "ERROR: -adda only works with -fixed orientation" << endl;
            exit(1);
        }
        additionalSummary += ", random orientation, ";
        AngleRange beta = GetRange(args, "b", particle);
        AngleRange gamma = GetRange(args, "g", particle);
        additionalSummary += "grid: " + to_string(beta.number) + "x" + to_string(gamma.number) + "\n\n";
        cout << additionalSummary;
        tracer.m_summary = additionalSummary;
        tracer.TraceRandom(beta, gamma);
    }

    delete handler;

    cout << endl << "done";

    if (!args.IsCatched("close"))
    {
        getchar();
    }

    return 0;
}
