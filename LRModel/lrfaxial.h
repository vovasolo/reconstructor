#ifndef LRFAXIAL_H
#define LRFAXIAL_H

#include "lrf.h"
#include <cmath>

class Bspline1d;
class BSfit1D;
class Compress1d;

class LRFaxial : public LRF
{
public:
    LRFaxial(double rmax, int nint);
    LRFaxial(const Json &json);
    LRFaxial(std::string &json_str);
    ~LRFaxial();

    virtual LRFaxial* clone() const;

    bool   inDomain(double x, double y, double z = 0) const override;
    double eval(double x, double y, double z = 0) const override;
    double getRmax() const override { return rmax; }

    double evalDrvX(double x, double y, double z = 0) const override;
    double evalDrvY(double x, double y, double z = 0) const override;

    bool fitData(const std::vector <LRFdata> &data) override;
    void addData(const std::vector <LRFdata> &data) override;
    bool doFit() override;

    std::string type() const override { return std::string("Axial"); }

    void ToJsonObject(Json_object &json) const override;

    // relative gain calculation
    double GetRatio(LRF* other) const override;

    //---

    virtual bool isReady () const;    

    int getNint() const { return nint; }
    //double fitRData(int npts, const double *r, const double *data);

    const Bspline1d *getSpline() const;

    //Andr: suggest all to lower case on start, add "get" to the lower ones
    void SetOrigin(double x0, double y0);
    void SetRmin(double rmin);
    void SetRmax(double rmax);
    void SetCompression(Compress1d *compress);

    void SetFlatTop(bool val) {flattop = val;}
    void SetNonIncreasing(bool val) {non_increasing = val;}

    // calculation of radius + provision for compression
    double R(double x, double y) const {return sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));}
    double R2(double x, double y) const {return (x-x0)*(x-x0)+(y-y0)*(y-y0);}
    double Rho(double r) const;
    double Rho(double x, double y) const;
    double RhoDrvX(double x, double y) const;
    double RhoDrvY(double x, double y) const;

protected:
    void Init();
    BSfit1D *InitFit();

protected:
    double x0   = 0;  // center
    double y0   = 0;  // center
    double rmin = 0;  // domain
    double rmax = 0;  // domain
    double rmin2;     // domain
    double rmax2;	  // domain
    int nint;		  // intervals

    bool flattop = false;   // set to true if you want to have zero derivative at the origin
    bool non_increasing = false;
    Bspline1d *bsr = nullptr; 	// spline describing radial dependence
    Compress1d *compress = nullptr; // optional compression
    BSfit1D *bsfit = nullptr;     // object used in fitting
    bool init_done = false;
    std::string json_err;
};

#endif // LRFAXIAL_H
