#include "reconstructor.h"
#include "lrmodel.h"
#include "TROOT.h"
#include <iostream>

static CostML *RecCostML = nullptr;
static CostChi2 *RecCostChi2 = nullptr;

Reconstructor::Reconstructor(LRModel *lrm)
{
    this->lrm = lrm;
    nsensors = lrm->GetSensorCount();
    sensor.resize(nsensors);
    active.resize(nsensors, true);
    for (int i=0; i<nsensors; i++) {
        sensor[i].x = lrm->GetX(i);
        sensor[i].y = lrm->GetY(i);
        sensor[i].gain = 1.;
        sensor[i].on = true;
    }
    A.resize(nsensors);
    sat.resize(nsensors);
}

bool Reconstructor::InitMinimizer()
{
    RootMinimizer = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
    //RootMinimizer = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kSimplex);
    RootMinimizer->SetMaxFunctionCalls(M2MaxFuncCalls);
    RootMinimizer->SetMaxIterations(M2MaxIterations);
    RootMinimizer->SetTolerance(M2Tolerance);

// Set Minuit2 and ROOT verbosity level
    RootMinimizer->SetPrintLevel(MinuitPrintLevel);
    gErrorIgnoreLevel = RootPrintLevel;

    if (method == ML) {
        RecCostML = new CostML(this);
        FunctorLSML = new ROOT::Math::Functor(*RecCostML, 3);
        LastMiniValue = 1.e100; // reset for the new event
    } else {
        RecCostChi2 = new CostChi2(this);
        FunctorLSML = new ROOT::Math::Functor(*RecCostChi2, 3);
        LastMiniValue = 1.e6; //reset for the new event
    }

    RootMinimizer->SetFunction(*FunctorLSML);
    return true;
}

// Check for static and dynamic passives + saturation
void Reconstructor::checkActive()
{
    nactive = 0;
    int maxid = getMaxSignalID();
    double cutoff = std::max(rec_abs_cutoff, A[maxid]*rec_rel_cutoff);
    for (int i=0; i<nsensors; i++) {
// sensor is active if it's enabled AND (NOT saturated) AND above the cutoff
        active[i] = sensor[i].on && !sat[i] && A[i] > cutoff;
// AND within cutoff radius
        active[i] = active[i] && getDistFromSensor(i, guess_x, guess_y) <= rec_cutoff_radius;
        nactive += active[i] ? 1 : 0;
    }
}

bool Reconstructor::ProcessEvent(std::vector <double> a, std::vector <bool> sat)
{
    if (a.size() < nsensors || sat.size() < nsensors)
        return false;
    for (int i=0; i<nsensors; i++) {
        A[i] = a[i]/sensor[i].gain;
        this->sat[i] = sat[i];
    }

// initial guess
    guessByCOG();
    guess_e = getSumSignal()*ecal;

// determine active sensors and see if there are enough for reconstruction
    checkActive();
    rec_dof = nactive - 3;
    if (rec_dof < 1) {
        rec_status = 6;
        return false;
    }

// set initial variables to minimize
    RootMinimizer->SetVariable(0, "x", guess_x, RMstepX);
    RootMinimizer->SetVariable(1, "y", guess_y, RMstepY);
    RootMinimizer->SetLowerLimitedVariable(2, "e", guess_e, guess_e*0.2, 1.0e-6);

    // do the minimization
    bool fOK = false;
    fOK = RootMinimizer->Minimize();

    if (fOK) {
        rec_status = 0 ;		// Reconstruction successfull

        const double *xs = RootMinimizer->X();
        rec_x = xs[0];
        rec_y = xs[1];
        rec_e = xs[2];
        rec_min = RootMinimizer->MinValue();

    // Calc Hessian matrix and get status
        int ndim = RootMinimizer->NDim();
        double cov[ndim*ndim];
        RootMinimizer->Hesse();
        RootMinimizer->GetCovMatrix(cov);
        cov_xx = cov[0]; // first column first row
        cov_yy = cov[ndim+1]; // second column second row
        cov_xy = cov[1];      // second column first row
        return true;
    } else {
        rec_status = RootMinimizer->Status(); // reason why it has failed
        return false;
    }
}

int Reconstructor::getMaxSignalID()
{
    auto strongest = std::max_element(std::begin(A), std::end(A));
    return std::distance(std::begin(A), strongest);
}

void Reconstructor::guessByMax()
{
    int maxid = getMaxSignalID();
    guess_x = sensor[maxid].x;
    guess_y = sensor[maxid].y;
}

void Reconstructor::guessByCOG()
{
    int maxid = getMaxSignalID();
    double cutoff = std::max(cog_abs_cutoff, A[maxid]*cog_rel_cutoff);
    double sum_x, sum_y, sum_dn;
    sum_x = sum_y = sum_dn = 0.;
    for (int i=0; i<nsensors; i++)
        if (A[i] >= cutoff) {
            sum_x += sensor[i].x*A[i];
            sum_y += sensor[i].y*A[i];
            sum_dn += A[i];
    }
    guess_x = sum_x/sum_dn;
    guess_y = sum_y/sum_dn;
}

double Reconstructor::getSumSignal()
{
    double sum = 0.;
    for (int i=0; i<nsensors; i++)
        sum += A[i];
    return sum;
}

double Reconstructor::getDistFromSensor(int id, double x, double y)
{
  return hypot(x-sensor[id].x, y-sensor[id].y);
}

double Reconstructor::getChi2(double x, double y, double z, double energy)
{
    double sum = 0;
    double r[3];
    r[0] = x; r[1] = y; r[2] = z;
//    std::cout << "X:" << x << "    Y:" << y << "   E:" << energy << std::endl;

    for (int i = 0; i < nsensors; i++) {
        if (!active[i])
            continue;

        double LRFhere = lrm->Eval(i, r)*energy; // LRF(X, Y, Z) * energy;
        if (LRFhere <= 0.)
            return LastMiniValue *= 1.25; //if LRFs are not defined for this coordinates

        double delta = (LRFhere - A[i]);
        sum += fWeightedLS ? delta*delta/LRFhere : delta*delta;
    }
//    std::cout << "Sum: " << sum << std::endl;
    return LastMiniValue = sum;
}

double Reconstructor::getLogLH(double x, double y, double z, double energy)
{
    double sum = 0;
    double r[3];
    r[0] = x; r[1] = y; r[2] = z;

    for (int i = 0; i < nsensors; i++) {
        if (!active[i])
            continue;

        double LRFhere = lrm->Eval(i, r)*energy; // LRF(X, Y, Z) * energy;
        if (LRFhere <= 0.)
            return LastMiniValue *= 1.25; //if LRFs are not defined for this coordinates

        sum += A[i]*log(LRFhere) - LRFhere; // measures probability
    }

    return LastMiniValue = sum;
}

double CostChi2::operator()(const double *p) // 0-x, 1-y, 2-energy
{
    return rec->getChi2(p[0], p[1], 0., p[2]);
}

double CostML::operator()(const double *p) // 0-x, 1-y, 2-energy
{
    return rec->getLogLH(p[0], p[1], 0., p[2]);
}
