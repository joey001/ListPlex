#pragma once
#include <stdlib.h>
#include <sys/time.h>
//gettimeofday
class Timer
{
private:
    double lastTime;
    //bool on;
    struct timezone tzp;
    double getTime()
    {
        timeval now;
        gettimeofday(&now, &tzp);
        return ((double)now.tv_sec) + ((double)now.tv_usec) / 1000000.;
    }

public:
    Timer()
    {
        struct timezone tz = {0, 0};
        tzp = tz;
        //on = false;
    }

    void start()
    {
        lastTime = getTime();
    }
    double stop()
    {
        double d = (getTime() - lastTime);
        return d;
    }
    void clear()
    {
        lastTime = 0.0;
    }
};