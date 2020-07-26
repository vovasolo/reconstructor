#ifndef COMPRESS_H
#define COMPRESS_H

#include <cmath>
#include "lrfio.h"

class Compress1d : public LRF_IO
{
public:
    Compress1d() {}
    virtual ~Compress1d() {}
    virtual Compress1d* clone() const = 0;
    virtual double Rho(double r) const = 0;
    virtual double RhoDrv(double r) const = 0;
    virtual void ToJsonObject(Json_object &json) const = 0;

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
    virtual DualSlopeCompress* clone() const { return new DualSlopeCompress(*this); }
    void Init();
    virtual double Rho(double r) const;
    virtual double RhoDrv(double r) const;
    virtual void ToJsonObject(Json_object &json) const;

private:
    double k;
    double r0;
    double lam;

    double a;
    double b;
    double lam2;
};

#endif // COMPRESS_H
