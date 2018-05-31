#pragma once

#include "string.H"
#include "clockTime.H"

namespace Foam
{
    class StopWatch
    {
        private:
            string _name;
            clockTime timer;
            double totalTime;
        public:
            StopWatch(string myname)
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
            string name() const
            {
                return this->_name;
            }
    };
}
