#ifndef STOPWATCH_H
#define STOPWATCH_H

#include "string.H"
#include "clockTime.H"

class Foam::StopWatch
{
    private:
        Foam::string _name;
        Foam::clockTime timer;
        double totalTime;
    public:
        Foam::StopWatch(string myname)
        {
            this->totalTime = 0;
            this->_name = myname;
        }
        void start()
        {
            // init new clock
            this->timer = Foam::clockTime();
        }
        void stop()
        {
            this->totalTime += this->timer.elapsedTime();
        }
        double getTotalTime() const
        {
            return this->totalTime;
        }
        Foam::string name() const
        {
            return this->_name;
        }
};

#endif
