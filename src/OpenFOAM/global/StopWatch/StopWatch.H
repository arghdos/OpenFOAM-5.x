#pragma once

#include <string.h>
#include "clockTime.H"

namespace Foam
{
    class StopWatch
    {
        private:
            std::string name;
            clockTime timer;
            double totalTime;
        public:
            StopWatch(std::string name)
            {
                this->totalTime = 0;
                this->name = name;
            }
            void start()
            {
                // init new clock
                this->timer = clockTime();
            }
            void stop()
            {
                this->totalTime += this->timer.elapsedTime();
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
}
