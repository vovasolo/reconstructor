#include "lrmodel.h"
#include "lrf.h"
#include "lrfaxial.h"
#include "transform.h"
#include "profileHist.h"
#include <cmath>
#include "json11.hpp"
//#include <complex>

double LRSensor::GetRadius() const
{
    return sqrt( x*x + y*y );
}

double LRSensor::GetPhi() const
{
    double phi = atan2 (y, x);
    // atan2 maps into [-pi;pi], make it [0;2*pi]
    if (phi < -1.0e-6) phi += 2.0*M_PI;
    return phi;
}

LRModel::LRModel(int n)
{
    Sensor.resize(n);
}

LRModel::LRModel(int n, LRF *default_lrf)
{
    Sensor.resize(n);
    DefaultLRF = default_lrf;
}

LRModel::~LRModel()
{
    ClearAll();
    delete DefaultLRF;
}

void LRModel::ClearAll()
{
    for (LRSensor s : Sensor) {
        delete s.tr;
        delete s.lrf;
    }
    for (LRGroup g : Group) {
        delete g.glrf;
    }
    Sensor.clear();
    Group.clear();
}

void LRModel::ResetGroups()
{
    for (int gid=0; gid<GetGroupCount(); gid++) {
        std::set <int> members = GroupMembers(gid);
        for (int i : members)
            RemoveFromGroup(i);
        delete Group.at(gid).glrf;
    }
    Group.clear();
}

void LRModel::AddSensor(int id, double x, double y)
{
    Sensor.at(id).id = id; // this is a canary to make sure that all sensors are initialized
    Sensor.at(id).x = x;
    Sensor.at(id).y = y;
//    Sensor.at(id).z = z;
    if (DefaultLRF)
        Sensor.at(id).lrf = DefaultLRF->clone();
}

std::vector <double> LRModel::GetAllX() const
{
    std::vector <double> tmp(GetSensorCount());
    for (LRSensor s : Sensor)
        tmp.push_back(s.x);
    return tmp;
}
std::vector <double> LRModel::GetAllY() const
{
    std::vector <double> tmp(GetSensorCount());
    for (LRSensor s : Sensor)
        tmp.push_back(s.y);
    return tmp;
}

int LRModel::CreateGroup()
{
    LRGroup g;
    if (DefaultLRF)
        g.glrf = DefaultLRF->clone();
    g.id = Group.size();
    Group.push_back(g);  
    return g.id;
}

bool LRModel::DissolveGroup(int gid)
{
    std::set <int> members = GroupMembers(gid);
    for (int i : members)
        RemoveFromGroup(i);

    delete Group.at(gid).glrf;
    Group.erase(Group.begin()+gid);
    return true;
}

bool LRModel::AddToGroup(int id, int gid, Transform *tr)
{
    if (!SensorExists(id) || !GroupExists(gid))
        return false;
    if (GetGroup(id) != -1)
        return false; // already belongs to another group

    SetGroup(id, gid);
    SetTransform(id, tr);
    GroupMembers(gid).insert(id);
    return true;
}

bool LRModel::RemoveFromGroup(int id, UngroupPolicy policy)
{
    int gid = GetGroup(id);
    if (gid == -1)
        return false; // no group => nothing to do

    SetGroup(id, -1);
    delete Sensor.at(id).lrf;
    switch (policy) {
        case KeepLRF:
            Sensor.at(id).lrf = Group.at(gid).glrf->clone();
            break;
        case ResetLRF:
            Sensor.at(id).lrf = DefaultLRF->clone();
            SetTransform(id, 0);
            SetGain(id, 1.0);
            break;
    }

    GroupMembers(gid).erase(id);
    if (GroupMembers(gid).size() == 0) {
        delete Group.at(gid).glrf;
        Group.erase(Group.begin()+gid);
    }
    return true;
}

void LRModel::SetTransform(int id, Transform *tr)
{
    delete Sensor.at(id).tr;
    Sensor.at(id).tr = tr;
}

void LRModel::SetLRF(int id, LRF *lrfptr)
{
    RemoveFromGroup(id);
    delete Sensor.at(id).lrf;
    Sensor.at(id).lrf = lrfptr;
}

LRF *LRModel::GetLRF(int id)
{
    int gid = GetGroup(id);
    return gid == -1 ? Sensor.at(id).lrf : Group.at(gid).glrf;
}

void LRModel::SetGroupLRF(int gid, LRF *lrfptr)
{
    delete Group.at(gid).glrf;
    Group.at(gid).glrf = lrfptr;
}

LRF *LRModel::GetGroupLRF(int gid)
{
    return Group.at(gid).glrf;
}

bool LRModel::InDomain(int id, double *pos_world)
{
    double x = pos_world[0];
    double y = pos_world[1];
    double z = pos_world[2];
    GetTransform(id)->DoTransform(&x, &y, &z);
    return GetLRF(id)->inDomain(x, y, z);
}

double LRModel::Eval(int id, double *pos_world)
{
    double x = pos_world[0];
    double y = pos_world[1];
    double z = pos_world[2];
    if (GetTransform(id))
        GetTransform(id)->DoTransform(&x, &y, &z);
    return GetLRF(id)->eval(x, y, z)*GetGain(id);
}

double LRModel::EvalDrvX(int id, double *pos_world)
{
    double x = pos_world[0];
    double y = pos_world[1];
    double z = pos_world[2];
    if (GetTransform(id))
        GetTransform(id)->DoTransform(&x, &y, &z);
    return GetLRF(id)->evalDrvX(x, y, z)*GetGain(id);
}

double LRModel::EvalDrvY(int id, double *pos_world)
{
    double x = pos_world[0];
    double y = pos_world[1];
    double z = pos_world[2];
    if (GetTransform(id))
        GetTransform(id)->DoTransform(&x, &y, &z);
    return GetLRF(id)->evalDrvY(x, y, z)*GetGain(id);
}

bool LRModel::FitNotBinnedData(int id, const std::vector <LRFdata> &data)
{
    std::vector <LRFdata> trdata;
    Transform *tr = GetTransform(id);
    double gain = GetGain(id);
    for (LRFdata d : data) {
        if (tr)
            tr->DoTransform(&d[0], &d[1], &d[2]);
        d[3] /= gain;
        trdata.push_back(d);
    }
    return GetLRF(id)->fitData(trdata);
}

void LRModel::AddFitData(int id, const std::vector <LRFdata> &data)
{
    std::vector <LRFdata> trdata;
    trdata.reserve(data.size());
    Transform *tr = GetTransform(id);
    double gain = GetGain(id);
    for (LRFdata d : data) {
        if (tr)
            tr->DoTransform(&d[0], &d[1], &d[2]);
        d[3] /= gain;
        trdata.push_back(d);
    }

    int gid = GetGroup(id);
    if (gid >= 0)
        GetGroupLRF(gid)->addData(trdata);
    else
        GetLRF(id)->addData(trdata);
}

bool LRModel::FitSensor(int id)
{
    return GetLRF(id)->doFit();
}

bool LRModel::FitGroup(int gid)
{
    return GetGroupLRF(gid)->doFit();
}

void LRModel::ClearAllFitData()
{

}

void LRModel::MakeGroupsCommon()
{
//    Reset();
    int gid = CreateGroup();
    for (LRSensor s : Sensor) {
        Transform *tr = new TranslateLRF(-s.x, -s.y);
        AddToGroup(s.id, gid, tr);
    }
}

void LRModel::MakeGroupsByRadius()
{
//    Reset();
    double R = 0.;
    std::vector <LRSensor> tmp_group;
    std::vector <LRSensor> local_sensors = Sensor;
    std::sort (local_sensors.begin(), local_sensors.end(), LRSensor::Compare_R);

    for (LRSensor s : local_sensors) {
        if (fabs(R-s.GetRadius()) <= tol)
            tmp_group.push_back(s);
        else {
            if (tmp_group.size() >= 2)
                MakeRotGroup(tmp_group);
            tmp_group.clear();
            tmp_group.push_back(s);
            R = s.GetRadius();
        }
    }

    if (tmp_group.size() >= 2)
        MakeRotGroup(tmp_group);
}

void LRModel::MakeRotGroup(std::vector <LRSensor> &ring)
{
    int gid = CreateGroup();
    std::sort (ring.begin(), ring.end(), LRSensor::Compare_Phi);
    double phi0 = ring[0].GetPhi();
    Group[gid].x = ring[0].x;
    Group[gid].y = ring[0].y;
    AddToGroup(ring[0].id, gid, 0);

    for (unsigned int i=1; i<ring.size(); i++)
        AddToGroup(ring[i].id, gid, new RotateLRF(phi0 - ring[i].GetPhi()));
}

void LRModel::MakeGroupsByTransform(std::vector <Transform*> vtr)
{
    // using local copy (ls) of the Sensor vector to simplify the housekeeping
    // shallow copy is OK as only the sensor positions are used
    std::vector <LRSensor> ls = Sensor;
    std::sort (ls.begin(), ls.end(), LRSensor::Compare_R);
    while (ls.size() > 1) {
        int gid = CreateGroup();
        AddToGroup(ls[0].id, gid, 0);
        double x0 = Group[gid].x = ls[0].x;
        double y0 = Group[gid].y = ls[0].y;
        for (Transform *tr : vtr) {
            for (unsigned int i=1; i<ls.size(); i++) {
                double x1 = ls[i].x;
                double y1 = ls[i].y;
                double z1 = 0.;
                tr->DoTransform(&x1, &y1, &z1); // must be forward transform
                if (fabs(x1-x0) < tol && fabs(y1-y0) < tol) {
                    AddToGroup(ls[i].id, gid, tr->clone());
                    ls.erase(ls.begin() + i);
                    break;
                }
            }
        }
        ls.erase(ls.begin());
        if (GetGroupMembersCount(gid) < 2)
            DissolveGroup(gid);
    }
}

std::vector <Transform*> LRModel::MakeVtrRectangle()
{
    std::vector <Transform*> vtr;
    vtr.push_back(new ReflectLRF(0));
    vtr.push_back(new ReflectLRF(M_PI/2));
    vtr.push_back(new RotateLRF(M_PI));
    return vtr;
}

std::vector <Transform*> LRModel::MakeVtrSquare()
{
    std::vector <Transform*> vtr;
    for (int i=0; i<4; i++)
        vtr.push_back(new ReflectLRF(i*M_PI/4));
    for (int i=1; i<4; i++)
        vtr.push_back(new RotateLRF(i*M_PI/2));
    return vtr;
}

std::vector <Transform*> LRModel::MakeVtrHexagon()
{
    std::vector <Transform*> vtr;
    for (int i=0; i<6; i++)
        vtr.push_back(new ReflectLRF(i*M_PI/6));
    for (int i=1; i<6; i++)
        vtr.push_back(new RotateLRF(i*M_PI/3));
    return vtr;
}

std::vector <Transform*> LRModel::MakeVtrNgon(int n)
{
    std::vector <Transform*> vtr;
    for (int i=0; i<n; i++)
        vtr.push_back(new ReflectLRF(i*M_PI/n));
    for (int i=1; i<n; i++)
        vtr.push_back(new RotateLRF(i*2*M_PI/n));
    return vtr;
}

// Input-output

Json::object LRModel::SensorGetJsonObject(int id) const
{
    Json::object json;
    LRSensor s = Sensor.at(id);
    json["id"] = s.id;
    json["group_id"] = s.group_id;
    json["x"] = s.x;
    json["y"] = s.y;
    json["gain"] = s.gain;
    if (s.tr) json["transform"] = s.tr->GetJsonObject();
    if (s.group_id != -1 && s.lrf) json["LRF"] = s.lrf->GetJsonObject();
    return json;
}

void LRModel::ReadSensor(const Json &json)
{
    int id = json["id"].int_value();
    LRSensor &s = Sensor.at(id);
    s.id = id;
    s.x = json["x"].number_value();
    s.y = json["y"].number_value();
    s.group_id = json["group_id"].int_value();
    s.gain = json["gain"].number_value();
    if (json["transform"].is_object())
        s.tr = Transform::Factory(json["transform"]);
    if (json["LRF"].is_object())
        s.lrf = new LRFaxial(json["LRF"]);
}

Json::object LRModel::GroupGetJsonObject(int gid) const
{
    Json::object json;
    LRGroup g = Group.at(gid);
    json["id"] = g.id;
    json["x0"] = g.x;
    json["y0"] = g.y;
    json["members"] = g.members;
    if (g.glrf) json["LRF"] = g.glrf->GetJsonObject();
    return json;
}

void LRModel::ReadGroup(const Json &json)
{
    int id = json["id"].int_value();
    LRGroup &g = Group.at(id);
    g.id = id;
    g.x = json["x"].number_value();
    g.y = json["y"].number_value();
    if (json["members"].is_array()) {
        Json::array members = json["members"].array_items();
        for (unsigned int i=0; i<members.size(); i++)
            g.members.insert(members[i].int_value());
    }
    if (json["LRF"].is_object())
        g.glrf = new LRFaxial(json["LRF"]);
}

void LRModel::ToJsonObject(Json_object &json) const
{
    std::vector <Json_object> sensors;
    std::vector <Json_object> groups;
    int n_sensors = (int)Sensor.size();
    int n_groups = (int)Group.size();
    json["n_sensors"] = n_sensors;
    json["n_groups"] = n_groups;

    for (int i=0; i<n_sensors; i++)
        sensors.push_back(SensorGetJsonObject(i));
    json["sensors"] = sensors;
    for (int i=0; i<n_groups; i++)
        groups.push_back(GroupGetJsonObject(i));
    json["groups"] = groups;
}

Json::object LRModel::GetJsonObject() const
{
    Json::object json;
    ToJsonObject(json);
    return json;
}

std::string LRModel::GetJsonString() const
{
    Json::object json;
    ToJsonObject(json);
    return Json(json).dump();
}

LRModel::LRModel(const Json &json)
{
    int n_sensors = json["n_sensors"].int_value();
    int n_groups = json["n_groups"].int_value();
    Sensor.resize(n_sensors);
    Group.resize(n_groups);
    if (json["sensors"].is_array()) {
        Json::array sensors = json["sensors"].array_items();
        for (unsigned int i=0; i<sensors.size(); i++)
            ReadSensor(sensors[i]);
    }
    if (json["groups"].is_array()) {
        Json::array groups = json["groups"].array_items();
        for (unsigned int i=0; i<groups.size(); i++)
            ReadGroup(groups[i]);
    }
}

LRModel::LRModel(std::string &json_str) : LRModel(Json::parse(json_str, json_err)) {}

// Utility
double LRModel::GetMaxR(int id, const std::vector <LRFdata> &data) const
{
    double x0 = GetX(id);
    double y0 = GetY(id);
    double maxr2 = 0.;

    for (LRFdata d : data) {
        double dx = d[0]-x0;
        double dy = d[1]-y0;
        maxr2 = std::max(maxr2, dx*dx + dy*dy);
    }
    return sqrt(maxr2);
}

double LRModel::GetGroupMaxR(int gid, const std::vector <LRFdata> &data) const
{
    double maxr = 0.;
    for (int id : Group.at(gid).members)
        maxr = std::max(maxr, GetMaxR(id, data));
    return maxr;
}

// Gain estimator

GainEstimator::GainEstimator(LRModel *M, int gid)
{
    this->M = M;
    this->gid = gid;
    members = M->GroupMembers(gid);
//    lrfs.reserve(members.size());
    for (int id : members) {
        lrfs.push_back(M->GetGroupLRF(gid)->clone());
    }
}

GainEstimator::~GainEstimator()
{
    for (int id : members)
        delete lrfs.at(id);
}

void GainEstimator::AddData(int id, const std::vector <LRFdata> &data)
{
    lrfs.at(id)->addData(data);
}

double GainEstimator::GetRelativeGain(int id, int refid)
{
    return lrfs.at(id)->GetRatio(lrfs.at(refid));
}

std::vector <double> GainEstimator::GetAllRelativeGains(int refid)
{
    std::vector <double> gains;
    for (unsigned int id = 0; id < members.size(); id++)
        gains.push_back(GetRelativeGain(id, refid));
    return gains;
}
