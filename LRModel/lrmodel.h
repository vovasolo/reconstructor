#ifndef LRMODEL_H
#define LRMODEL_H

#include <vector>
#include <set>
#include <cmath>
#include "lrf.h"
#include "lrfio.h"
#include "profileHist.h"

class Transform;

class LRSensor
{
friend class LRModel;
public:
    double GetRadius() const;
    double GetPhi() const;
    double GetDistance(double x1, double y1) const
        {return sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));}
    bool operator < (const LRSensor& s) const {
        return (GetRadius() < s.GetRadius());
    }
    static bool Compare_R(const LRSensor &a, const LRSensor &b)
        { return (a.GetRadius() < b.GetRadius()); }
    static bool Compare_Phi(const LRSensor &a, const LRSensor &b)
        { return (a.GetPhi() < b.GetPhi()); }
    static double Distance(const LRSensor &a, const LRSensor &b)
        { return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)); }

protected:
    int id = -1;        // must coincide with vector element
    int group_id = -1;  // no group by default
    double x, y;        // PMT position
//  double z = 0.;      // later
//  double normal[3];   // later
    double gain = 1.0;  // relative gain
    Transform *tr = 0;
    LRF *lrf = 0;
};

struct LRGroup
{
    int id = -1;
    std::set <int> members;
    double x, y;        // coordinates of the reference point
//  double z = 0.;      // later
//  double normal[3];   // later
    LRF *glrf = 0;
};

class LRModel
{
    enum UngroupPolicy {
        KeepLRF,       // clone the group LRF and keep it
        ResetLRF       // clone the default LRF
    };

public:
    LRModel(int n);
    LRModel(int n, LRF *default_lrf);
    ~LRModel();

    void SetGain(int id, double gain) {Sensor.at(id).gain = gain;}
    void SetGroup(int id, int group) {Sensor.at(id).group_id = group;}
//    void SetLRFid(int id, int lrfid) {Sensor.at(id).lrfid = lrfid;}
//    void SetGroupLRFid(int gid, int lrfid) {Group.at(gid).lrfid = lrfid;}
    void SetTransform(int id, Transform *tr);

    double GetGain(int id) const {return Sensor.at(id).gain;}
    int GetGroup(int id) const {return Sensor.at(id).group_id;}
//    int  GetLRFid(int id) const {return Sensor.at(id).lrfid;}
//    int  GetGroupLRFid(int gid) const {return Group.at(gid).lrfid;}
    Transform *GetTransform(int id) const {return Sensor.at(id).tr;}
    double GetX(int id) const {return Sensor.at(id).x;}
    double GetY(int id) const {return Sensor.at(id).y;}
    std::vector <double> GetAllX() const;
    std::vector <double> GetAllY() const;
    double GetGroupX(int gid) const {return Group.at(gid).x;}
    double GetGroupY(int gid) const {return Group.at(gid).y;}
    double GetRadius(int id) const {return Sensor.at(id).GetRadius();}
    double GetPhi(int id) const {return Sensor.at(id).GetPhi();}
    double GetDistance(int ida, int idb) const
        {return LRSensor::Distance(Sensor.at(ida), Sensor.at(idb));}

    std::set <int> &GroupMembers(int gid) {return Group.at(gid).members;}

    int GetSensorCount() const {return Sensor.size();}
    int GetGroupCount() const {return Group.size();}
    int GetGroupMembersCount(int gid) const {return Group.at(gid).members.size();}
    bool SensorExists(int id) const {return id>=0 && id<Sensor.size();}
    bool GroupExists(int gid) const {return gid>=0 && gid<Group.size();}
//    int GetLRFCount() const {return Lrf.size();}

    void ClearAll();
    void ResetGroups();
    void AddSensor(int id, double x, double y);

    int CreateGroup();
    bool DissolveGroup(int gid);
    bool AddToGroup(int id, int group, Transform *tr);
    bool RemoveFromGroup(int id, UngroupPolicy policy = KeepLRF);

    void MakeGroupsCommon();
    void MakeGroupsByRadius();
    void MakeGroupsRectangle() { MakeGroupsByTransform(MakeVtrRectangle()); }
    void MakeGroupsSquare() { MakeGroupsByTransform(MakeVtrSquare()); }
    void MakeGroupsHexagon() { MakeGroupsByTransform(MakeVtrHexagon()); }
    void MakeGroupsNgon(int n) { MakeGroupsByTransform(MakeVtrNgon(n)); }

    void MakeRotGroup(std::vector <LRSensor> &ring);
    void MakeGroupsByTransform(std::vector <Transform*> vtr);
    std::vector <Transform*> MakeVtrRectangle();
    std::vector <Transform*> MakeVtrSquare();
    std::vector <Transform*> MakeVtrHexagon();
    std::vector <Transform*> MakeVtrNgon(int n);
    void SetTolerance(double tolerance) {tol = tolerance;}

// Access to LRFs
    void SetLRF(int id, LRF *lrfptr);
    LRF *GetLRF(int id);
    void SetGroupLRF(int gid, LRF *lrfptr);
    LRF *GetGroupLRF(int gid);
    void SetDefaultLRF(LRF *default_lrf) {DefaultLRF = default_lrf;}

// Evaluation
    bool InDomain(int id, double *pos_world);
    double Eval(int id, double *pos_world);
    double EvalLocal(int id, double *pos_local) { return GetLRF(id)->eval(pos_local)*GetGain(id); }
    double EvalDrvX(int id, double *pos_world);
    double EvalDrvY(int id, double *pos_world);

// Fitting
    // direct (not binned)
    bool FitNotBinnedData(int id, const std::vector <LRFdata> &data);
    // binned
    void AddFitData(int id, const std::vector <LRFdata> &data);
    bool FitSensor(int id);
    bool FitGroup(int gid);
    void ClearAllFitData();

// Save and Load
    Json_object SensorGetJsonObject(int id) const;
    Json_object GroupGetJsonObject(int gid) const;
    void ToJsonObject(Json_object &json) const;
    Json_object GetJsonObject() const;
    std::string GetJsonString() const;

    void ReadSensor(const Json &json);
    void ReadGroup(const Json &json);
    LRModel(const Json &json);
    LRModel(std::string &json_str);

// Utility
    double GetMaxR(int id, const std::vector <LRFdata> &data) const;
    double GetGroupMaxR(int gid, const std::vector <LRFdata> &data) const;

protected:
    std::vector <LRSensor> Sensor;
    std::vector <LRGroup> Group;
//    std::vector <LRF*> Lrf;
    LRF *DefaultLRF = 0;
    std::string json_err;
    double tol = 1.0e-4;
};

class GainEstimator
{
public:
    GainEstimator(LRModel *M, int gid);
    ~GainEstimator();
    void AddData(int id, const std::vector <LRFdata> &data);
    double GetRelativeGain(int id, int refid);
    std::vector <double> GetAllRelativeGains(int refid);

protected:
    LRModel *M;
    int gid;
    std::vector <LRF*> lrfs;
    std::set <int> members;
};

#endif // LRMODEL_H
