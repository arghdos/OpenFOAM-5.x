//declare timers
#include "StopWatch.H"
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
//- the internal clock for qdot eval
extern StopWatch QDotEvalTime;
extern StopWatch SetYInertTime;
extern StopWatch YLoopTime;
//- chemical time scale eval
extern StopWatch TCEvalTime;

// ODE-solver timers
extern StopWatch ODESolveTime;
extern StopWatch JacobianEvalTime;
extern StopWatch dYdTEvalTime;
