TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += lib
INCLUDEPATH += spline123
INCLUDEPATH += LRModel
INCLUDEPATH += /usr/include/eigen3
INCLUDEPATH += $$system(root-config --incdir)

LIBS += $$system(root-config --libs) -lGeom -lGeomPainter -lGeomBuilder -lMinuit2 -lSpectrum -ltbb

SOURCES += \
    LRModel/lrmodel.cpp \
    LRModel/lrfaxial.cpp \
    LRModel/transform.cpp \
    LRModel/compress.cpp \
    LRModel/lrfio.cpp \
    spline123/bsfit123.cpp \
    spline123/profileHist.cpp \
    spline123/bspline123d.cpp \
    lib/json11.cpp \
    reconstructor.cpp \
    example1.cpp

HEADERS += \
    LRModel/lrfaxial.h \
    LRModel/lrmodel.h \
    LRModel/compress.h \
    LRModel/lrf.h \
    LRModel/transform.h \
    LRModel/lrfio.h \
    spline123/bsfit123.h \
    spline123/profileHist.h \
    spline123/bspline123d.h \
    lib/json11.hpp \
    lib/eiquadprog.hpp \
    reconstructor.h
