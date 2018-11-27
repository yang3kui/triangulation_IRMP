#pragma once
#include <inttypes.h>
#include <chrono>
namespace tri
{
using time_t = uint64_t;
// clock in microseconds
class Clock
{
public:
    void reset();
    time_t now();
private:


    std::chrono::steady_clock::time_point start_;

};
}
