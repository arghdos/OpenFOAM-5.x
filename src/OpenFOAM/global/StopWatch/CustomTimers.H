#ifndef CUSTOMTIMERS_H
#define CUSTOMTIMERS_H

#include "StopWatch.H"

//declare timers
namespace Foam
{
    extern StopWatch totalTime;
    extern StopWatch mainLoopTime;
    extern StopWatch readControlsTime;
    extern StopWatch setDeltaTTime;
    extern StopWatch pimpleTime;
    extern StopWatch RhoEqnTime;
    extern StopWatch UEqnTime;
    extern StopWatch YEqnTime;
    extern StopWatch EEqnTime;
    extern StopWatch pEqnTime;
    extern StopWatch turbEqnTime;
    extern StopWatch writeTime;
    extern StopWatch infoTime;
    extern StopWatch rhoFetchTime;
    // equation specific timers
    // Y-eqn
    extern StopWatch YConvectionTime;
    extern StopWatch CombustionModelTime;
    extern StopWatch HeatReleaseTime;
    extern StopWatch SetYInertTime;
    extern StopWatch YLoopTime;
    //- chemistry model timers
    //- the internal clock for qdot eval
    extern StopWatch QDotEvalTime;
    //- chemical time scale eval
    extern StopWatch TCEvalTime;
    // ODE-solver timers
    extern StopWatch ODESolveTime;
    extern StopWatch JacobianEvalTime;
    extern StopWatch dYdTEvalTime;
    extern StopWatch ReactionRateEvalTime;
    extern StopWatch OmegaEvalTime;
    // combustion model timers
    extern StopWatch CombustionCorrectTime;
    extern StopWatch CombustionHeatReleaseTime;
    extern StopWatch CombustionFuelConsumptionTime;
}

#endif
