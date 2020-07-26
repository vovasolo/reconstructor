#include "compress.h"
#include "json11.hpp"

Compress1d* Compress1d::Factory(const Json &json)
{
    if (json["method"].is_string() && json["method"].string_value() == "dualslope") {
        Compress1d *comp = new DualSlopeCompress(json);
        if (comp->fValid)
            return comp;
        else {
            delete comp;
            return NULL;
        }
    }
    else
        return NULL;
}

void DualSlopeCompress::Init()
{
    if (r0 < 0. || k <= 1.)
        return;

    a = (k+1)/(k-1);
    lam2 = lam*lam;
    b = sqrt(r0*r0+lam2)+a*r0;

    fValid = true;
}

DualSlopeCompress::DualSlopeCompress(double k, double r0, double lam)
{
    this->r0 = r0;
    this->k = k;
    this->lam = lam;

    Init();
}

DualSlopeCompress::DualSlopeCompress(const Json &json)
{
    if (!json["r0"].is_number() || !json["k"].is_number() || !json["lam"].is_number())
        return;

    r0 = json["r0"].number_value();
    k = json["k"].number_value();
    lam = json["lam"].number_value();

    Init();
}

double DualSlopeCompress::Rho(double r) const
{
    double dr = r - r0;
    return std::max(0., b + dr*a - sqrt(dr*dr + lam2));
}

double DualSlopeCompress::RhoDrv(double r) const
{
    double dr = r - r0;
    return dr/sqrt(dr*dr + lam2) + a;
}

void  DualSlopeCompress::ToJsonObject(Json_object &json) const
{
    json["method"] = "dualslope";
    json["r0"] = r0;
    json["lam"] = lam;
    json["k"] = k;
}

