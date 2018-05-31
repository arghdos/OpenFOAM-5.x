#pragma once

#include <ctime>
#include <string.h>

class StopWatch
{
    private:
        std::clock_t startTime;
        std::string name;
        double totalTime;
    public:
        StopWatch(std::string name)
        {
            this->totalTime = 0;
            this->name = name;
        }
        void start()
        {
            this->startTime = std::clock();
        }
        void stop()
        {
            this->totalTime += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
        }
        double getTotalTime() const
        {
            return this->totalTime;
        }
        std::str name() const
        {
            return this->name;
        }
};
