//declare timers
#include "StopWatch.H"
extern Foam::StopWatch totalTime;
extern Foam::StopWatch mainLoopTime;
extern Foam::StopWatch readControlsTime;
extern Foam::StopWatch setDeltaTTime;
extern Foam::StopWatch pimpleTime;
extern Foam::StopWatch RhoEqnTime;
extern Foam::StopWatch UEqnTime;
extern Foam::StopWatch YEqnTime;
extern Foam::StopWatch EEqnTime;
extern Foam::StopWatch pEqnTime;
extern Foam::StopWatch turbEqnTime;
extern Foam::StopWatch writeTime;
extern Foam::StopWatch infoTime;
extern Foam::StopWatch rhoFetchTime;
// equation specific timers
// Y-eqn
extern Foam::StopWatch YConvectionTime;
extern Foam::StopWatch CombustionModelTime;
extern Foam::StopWatch HeatReleaseTime;
//- the internal clock for qdot eval
extern Foam::StopWatch QDotEvalTime;
extern Foam::StopWatch SetYInertTime;
extern Foam::StopWatch YLoopTime;
//- chemical time scale eval
extern Foam::StopWatch TCEvalTime;

// ODE-solver timers
extern Foam::StopWatch ODESolveTime;
extern Foam::StopWatch JacobianEvalTime;
extern Foam::StopWatch dYdTEvalTime;
