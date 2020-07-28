#ifndef COMPRESS_H
#define COMPRESS_H

#include <cmath>
#include "lrfio.h"

class Compress1d : public LRF_IO
{
public:
    //Compress1d() {}
    //virtual ~Compress1d() {} //Andr: already virtual in the base class

    virtual Compress1d* clone() const = 0;
    virtual double Rho(double r) const = 0;
    virtual double RhoDrv(double r) const = 0;

    void ToJsonObject(Json_object &json) const override = 0;

    static Compress1d* Factory(const Json &json);

protected:
    bool fValid = false;
};

// Dual slope

class DualSlopeCompress : public Compress1d
{
public:
    DualSlopeCompress(double k, double r0, double lam);
    DualSlopeCompress(const Json &json);
    DualSlopeCompress* clone() const override { return new DualSlopeCompress(*this); }
    double Rho(double r) const override;
    double RhoDrv(double r) const override;
    void ToJsonObject(Json_object &json) const override;

protected:
    void Init();

private:
    double k;
    double r0;
    double lam = 0;

    double a;
    double b;
    double lam2;
};

#endif // COMPRESS_H
