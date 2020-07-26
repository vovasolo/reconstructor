#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <vector>
#include "TMath.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

class LRModel;

struct RecSensor
{
    double x;
    double y;
    double gain;
    bool on;
};

class Reconstructor
{
public:
    enum Method {
        LS,            // least squares
        ML,            // maximum likelyhood
    };

public:
    Reconstructor(LRModel *lrm);
    ~Reconstructor() {;}

    bool InitMinimizer();
    bool ProcessEvent(std::vector <double> a, std::vector <bool> sat);

// cost functions
    double getChi2(double x, double y, double z, double energy);
    double getLogLH(double x, double y, double z, double energy);

// tracking of minimized value
    double LastMiniValue;

// public interface
public:
    double getGuessX() {return guess_x;}
    double getGuessY() {return guess_y;}
    double getGuessE() {return guess_e;}
    int getRecStatus() {return rec_status;}
    int getDof() {return rec_dof;}
    double getRecX() {return rec_x;}
    double getRecY() {return rec_y;}
    double getRecE() {return rec_e;}
    double getRecMin() {return rec_min;}
    double getChi2Min() {return getChi2(rec_x, rec_y, 0., rec_e);}
    double getCovXX() {return cov_xx;}
    double getCovYY() {return cov_yy;}
    double getCovXY() {return cov_xy;}
    void setCogAbsCutoff(double val) {cog_abs_cutoff = val;}
    void setCogRelCutoff(double val) {cog_rel_cutoff = val;}
    void setRecAbsCutoff(double val) {rec_abs_cutoff = val;}
    void setRecRelCutoff(double val) {rec_rel_cutoff = val;}
    void setRecCutoffRadius(double val) {rec_cutoff_radius = val;}
    void setEnergyCalibration(double val) {ecal = val;}
    void setGain(int id, double val) {sensor.at(id).gain = val;}

protected:
    LRModel *getLRModel() {return lrm;}
    void checkActive();
    double getSumSignal();
    int getMaxSignalID();
    void guessByMax();
    void guessByCOG();
    double getDistFromSensor(int id, double x, double y);

protected:
    LRModel *lrm;
    int nsensors = 0;
    int nactive = 0;
// cached sensor parameters
    std::vector <RecSensor> sensor;
    std::vector <bool> active;
// cached input parameters
    std::vector <double> A;
    std::vector <bool> sat;

// CoG
    double cog_abs_cutoff = 0.;
    double cog_rel_cutoff = 0.;

// initial guess
    double guess_x;
    double guess_y;
    double guess_e;
    double ecal = 3.75e-5; // approximate scaling factor between SumSignal and energy

// dynamic passives
    double rec_cutoff_radius = 1.0e12; // all by default
    double rec_abs_cutoff = 0.;
    double rec_rel_cutoff = 0.;

// method
    Method method = LS;
    bool fWeightedLS = true;

// ROOT/Minuit stuff
    ROOT::Math::Functor *FunctorLSML;
    ROOT::Minuit2::Minuit2Minimizer *RootMinimizer;
// initial steps
    double RMstepX = 1.;
    double RMstepY = 1.;
    double RMstepEnergy;
// control over MINUIT2 stopping conditions
    int M2MaxFuncCalls = 500;       // Max function calls
    int M2MaxIterations = 1000; 	// Max iterations
    double M2Tolerance = 0.001;		// Iteration stops when the function is within <tolerance> from the (estimated) min/max
// control over ROOT/MINUIT2 output
    int MinuitPrintLevel = 0;       // MINUIT2 messages
    int RootPrintLevel = 1001;      // ROOT messsages

// reconstruction result
    int rec_status;         // returned status of reconstruction
    double rec_x;			// reconstructed X position
    double rec_y;			// reconstructed Y position
    double rec_e;           // reconstructed energy
    double rec_min;         // reduced best chi-squared from reconstruction
    int rec_dof;			// degrees of freedom
    double cov_xx;		// variance in x
    double cov_yy;		// variance in y
    double cov_xy;		// covariance xy
};

class CostChi2
{
    public:
        CostChi2(Reconstructor *r) : rec(r) {;}
        double operator()(const double *p);
    private:
        Reconstructor *rec;
};

class CostML
{
    public:
        CostML(Reconstructor *r) : rec(r) {;}
        double operator()(const double *p);
    private:
        Reconstructor *rec;
};

#endif // RECONSTRUCTOR_H
