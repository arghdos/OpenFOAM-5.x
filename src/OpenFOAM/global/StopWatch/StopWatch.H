#pragma once

#include <string.h>
#include "clockTime.H"

namespace Foam
{
    class StopWatch
    {
        private:
            std::string _name;
            clockTime timer;
            double totalTime;
        public:
            StopWatch(std::string myname)
            {
                this->totalTime = 0;
                this->_name = myname;
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
            std::string name() const
            {
                return this->_name;
            }
    };
}
