#include "clock.h"
namespace tri
{
void Clock::reset()
{
    start_ = std::chrono::steady_clock::now();
}

time_t Clock::now()
{
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
}
}
