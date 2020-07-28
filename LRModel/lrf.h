#ifndef LRF_H
#define LRF_H

#include "lrfio.h"
#include <array>
#include <vector>
#include <string>

typedef std::array <double, 4> LRFdata;

class BSfit;

class LRF : public LRF_IO
{
public:
    //LRF() {}
    //virtual ~LRF() {} // Already virtual in the base class

    virtual LRF* clone() const = 0;

    bool   inDomain(double *pos) const {return inDomain(pos[0], pos[1], pos[2]);} // Andr: propose to change to   isInDomain()
    double eval(double *pos)     const {return eval(pos[0], pos[1], pos[2]);}
    double evalDrvX(double *pos) const {return evalDrvX(pos[0], pos[1], pos[2]);}
    double evalDrvY(double *pos) const {return evalDrvY(pos[0], pos[1], pos[2]);}

    void   setBinned(bool val) {binned = val;} // !*! Andr: invalidate ready flags?
    void   setNonNegative(bool val) {non_negative = val;}

    virtual bool inDomain(double x, double y, double z = 0) const = 0;  // Andr: propose to change to   isInDomain()
    // propose: all user interface to make non-virtual, with protected virtual functions inside - the same as inDomain(double *pos)

    virtual double eval(double x, double y, double z = 0) const = 0;
    virtual double evalDrvX(double x, double y, double z = 0) const = 0;
    virtual double evalDrvY(double x, double y, double z = 0) const = 0;

    virtual bool fitData(const std::vector <LRFdata> &data) = 0;
    virtual void addData(const std::vector <LRFdata> &data) = 0;
    virtual bool doFit() = 0;

    virtual std::string type() const = 0;  // Andr: rename to getType
    virtual bool isValid() const { return valid; }
//    virtual bool isReady () const;

    virtual double getRmax() const = 0;
    virtual double getXmin() const {return xmin;}
    virtual double getXmax() const {return xmax;}
    virtual double getYmin() const {return ymin;}
    virtual double getYmax() const {return ymax;}
    virtual double getZmin() const {return zmin;}
    virtual double getZmax() const {return zmax;}

    virtual double GetRatio(LRF* other) const = 0;   // Andr: suggest Get -> get...

protected:
    bool valid = false; // indicates if the LRF can be used for reconstruction
    double xmin, xmax; 	// xrange
    double ymin, ymax; 	// yrange
    double zmin, zmax; 	// zrange
    bool binned = true;
    bool non_negative = false;
};

#endif // LRF_H
