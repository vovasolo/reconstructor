#ifndef LRFIO_H
#define LRFIO_H

// to remove serialization routines and dependence on json11 library
// comment the following line
#define BSIO

#include <map>
#include <string>

namespace json11 {
    class Json;
}
using json11::Json;
typedef std::map<std::string, Json> Json_object;

class LRF_IO
{
public:
    //LRF_IO() {}
    virtual ~LRF_IO() {}

    virtual void ToJsonObject(Json_object &json) const = 0;   // Andr: suggest To... -> to...

    Json_object GetJsonObject() const; // get...
    std::string GetJsonString() const; // get...
};

#endif // LRFIO_H
