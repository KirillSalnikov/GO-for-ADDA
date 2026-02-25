CXX       = g++
CXXFLAGS  = -pipe -O2 -std=c++11 -msse4.2 -Wall -Wextra -fPIC
LDFLAGS   = -Wl,-O1
LIBS      = -lpthread -lm

SRCDIR    = src
BINDIR    = bin
TARGET    = $(BINDIR)/mbs

INCPATH   = -I$(SRCDIR) \
            -I$(SRCDIR)/math \
            -I$(SRCDIR)/bigint \
            -I$(SRCDIR)/handler \
            -I$(SRCDIR)/service \
            -I$(SRCDIR)/geometry \
            -I$(SRCDIR)/geometry/sse \
            -I$(SRCDIR)/splitting \
            -I$(SRCDIR)/common \
            -I$(SRCDIR)/particle \
            -I$(SRCDIR)/scattering \
            -I$(SRCDIR)/tracer \
            -I$(SRCDIR)/adda

SOURCES   = $(SRCDIR)/Beam.cpp \
            $(SRCDIR)/CalcTimer.cpp \
            $(SRCDIR)/main.cpp \
            $(SRCDIR)/ScatteringFiles.cpp \
            $(SRCDIR)/Splitting.cpp \
            $(SRCDIR)/Tracer.cpp \
            $(SRCDIR)/Tracks.cpp \
            $(SRCDIR)/adda/ADDAField.cpp \
            $(SRCDIR)/common/common.cpp \
            $(SRCDIR)/common/Matrix4x4.cpp \
            $(SRCDIR)/common/MullerMatrix.cpp \
            $(SRCDIR)/geometry/Facet.cpp \
            $(SRCDIR)/geometry/geometry_lib.cpp \
            $(SRCDIR)/geometry/Intersection.cpp \
            $(SRCDIR)/geometry/Polygon.cpp \
            $(SRCDIR)/handler/Handler.cpp \
            $(SRCDIR)/handler/HandlerGO.cpp \
            $(SRCDIR)/handler/HandlerTotalGO.cpp \
            $(SRCDIR)/handler/HandlerTracksGO.cpp \
            $(SRCDIR)/math/compl.cpp \
            $(SRCDIR)/math/JonesMatrix.cpp \
            $(SRCDIR)/math/matrix.cpp \
            $(SRCDIR)/math/Mueller.cpp \
            $(SRCDIR)/math/PhysMtr.cpp \
            $(SRCDIR)/particle/Bullet.cpp \
            $(SRCDIR)/particle/BulletRosette.cpp \
            $(SRCDIR)/particle/CertainAggregate.cpp \
            $(SRCDIR)/particle/ConcaveHexagonal.cpp \
            $(SRCDIR)/particle/Droxtal.cpp \
            $(SRCDIR)/particle/Hexagonal.cpp \
            $(SRCDIR)/particle/HexagonalAggregate.cpp \
            $(SRCDIR)/particle/Particle.cpp \
            $(SRCDIR)/particle/TiltedHexagonal.cpp \
            $(SRCDIR)/scattering/Scattering.cpp \
            $(SRCDIR)/scattering/ScatteringConvex.cpp \
            $(SRCDIR)/scattering/ScatteringNonConvex.cpp \
            $(SRCDIR)/tracer/TracerGO.cpp \
            $(SRCDIR)/geometry/intrinsic/intrinsics.cpp \
            $(SRCDIR)/bigint/BigInteger.cc \
            $(SRCDIR)/bigint/BigIntegerAlgorithms.cc \
            $(SRCDIR)/bigint/BigIntegerUtils.cc \
            $(SRCDIR)/bigint/BigUnsigned.cc \
            $(SRCDIR)/bigint/BigUnsignedInABase.cc

OBJECTS   = $(patsubst $(SRCDIR)/%.cpp,build/%.o,$(filter %.cpp,$(SOURCES))) \
            $(patsubst $(SRCDIR)/%.cc,build/%.o,$(filter %.cc,$(SOURCES)))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(BINDIR):
	mkdir -p $(BINDIR)

build/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

build/%.o: $(SRCDIR)/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

clean:
	rm -rf build $(TARGET)
