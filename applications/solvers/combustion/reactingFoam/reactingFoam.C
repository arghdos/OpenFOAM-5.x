/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "CustomTimers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void print_timers(List<Foam::StopWatch>& const watches, std::bool normalize=false)
{
    totalTime = watches.begin().getTotalTime();
    Info<<"Time Profile: ";
    for (List<Foam::StopWatch>::const_iterator it = watches.begin() + 1; it != watches.end(); ++it)
    {
        double time = it->getTotalTime();
        if (normalize)
        {
            time = (time / totalTime) * 100.0;
        }
        Info<<"\n\t" << it->name << " (s):" << time;
    }
    Info<<endl;
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //declare timers
    Foam::StopWatch totalTime(Foam::string("total runtime"));
    Foam::StopWatch mainLoopTime(Foam::string("main time-loop"));
    Foam::StopWatch readControlsTime(Foam::string("read controls"));
    Foam::StopWatch setDeltaTTime(Foam::string("set timestep"));
    Foam::StopWatch pimpleTime(Foam::string("pimple loop"));
    Foam::StopWatch RhoEqnTime(Foam::string("density equation"));
    Foam::StopWatch UEqnTime(Foam::string("velocity equations"));
    Foam::StopWatch YEqnTime(Foam::string("species equations"));
    Foam::StopWatch EEqnTime(Foam::string("energy equation"));
    Foam::StopWatch pEqnTime(Foam::string("pressure equation"));
    Foam::StopWatch turbEqnTime(Foam::string("turbulence equation"));
    Foam::StopWatch writeTime(Foam::string("write time"));
    Foam::StopWatch rhoFetchTime(Foam::string("read density"));
    Foam::StopWatch YConvectionTime(Foam::string("species convection initialization"));
    Foam::StopWatch CombustionModelTime(Foam::string("combustion mode evaluation"));
    Foam::StopWatch HeatReleaseTime(Foam::string("(outer) heat release evaluation"));
    Foam::StopWatch SetYInertTime(Foam::string("inert species handling"));
    Foam::StopWatch YLoopTime(Foam::string("species loop solution"));
    List<Foam::StopWatch> TimerList({
        totalTime,
        mainLoopTime,
        readControlsTime,
        setDeltaTTime,
        pimpleTime,
        RhoEqnTime,
        UEqnTime,
        YEqnTime,
        EEqnTime,
        pEqnTime,
        turbEqnTime,
        writeTime,
        rhoFetchTime,
        YConvectionTime,
        CombustionModelTime,
        HeatReleaseTime,
        YLoopTime,
        TCEvalTime,
        ODESolveTime,
        JacobianEvalTime,
        dYdTEvalTime,
        ReactionRateEvalTime,
        OmegaEvalTime
    });
    //end

    totalTime.start();

    while (runTime.run())
    {
        mainLoopTime.start();

        readControlsTime.start();
        #include "readTimeControls.H"
        readControlsTime.stop();

        setDeltaTTime.start()
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }
        setDeltaTTime.stop();

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        RhoEqnTime.start();
        #include "rhoEqn.H"
        RhoEqnTime.stop();

        pimpleTime.start();
        while (pimple.loop())
        {
            UEqnTime.start();
            #include "UEqn.H"
            UEqnTime.stop();

            YEqnTime.start();
            #include "YEqn.H"
            YEqnTime.stop();
            EEqnTime.start();
            #include "EEqn.H"
            EEqnTime.stop();

            // --- Pressure corrector loop
            pEqnTime.start();
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
            }
            pEqnTime.stop();

            turbEqnTime.start();
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            turbEqnTime.stop();
        }
        pimpleTime.stop();

        rhoFetchTime.start();
        rho = thermo.rho();
        rhoFetchTime.end();

        writeTime.start();
        runTime.write();
        writeTime.stop();

        print_timers(TimerList, false);
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    totalTime.stop();
    print_timers(TimerList, true);
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
