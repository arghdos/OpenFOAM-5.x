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

void print_timers(List<const StopWatch*> &watches, bool normalize=false)
{
    double totalTime = (*watches.begin())->getTotalTime();
    Info<<"Time Profile: ";
    for (List<const StopWatch*>::const_iterator it = watches.begin() + 1; it != watches.end(); ++it)
    {
        double time = (*it)->getTotalTime();
        if (normalize)
        {
            time = (time / totalTime) * 100.0;
            Info<<"\n\t" << (*it)->name() << " (%):" << time;
        }
        else
        {
            Info<<"\n\t" << (*it)->name() << " (s):" << time;
        }
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
    StopWatch totalTime(string("total runtime"));
    StopWatch mainLoopTime(string("main time-loop"));
    StopWatch readControlsTime(string("read controls"));
    StopWatch setDeltaTTime(string("set timestep"));
    StopWatch pimpleTime(string("pimple loop"));
    StopWatch RhoEqnTime(string("density equation"));
    StopWatch UEqnTime(string("velocity equations"));
    StopWatch YEqnTime(string("species equations"));
    StopWatch EEqnTime(string("energy equation"));
    StopWatch pEqnTime(string("pressure equation"));
    StopWatch turbEqnTime(string("turbulence equation"));
    StopWatch writeTime(string("write time"));
    StopWatch rhoFetchTime(string("read density"));
    StopWatch YConvectionTime(string("species convection initialization"));
    StopWatch CombustionModelTime(string("combustion mode evaluation"));
    StopWatch HeatReleaseTime(string("(outer) heat release evaluation"));
    StopWatch SetYInertTime(string("inert species handling"));
    StopWatch YLoopTime(string("species loop solution"));
    List<const StopWatch*> TimerList({
        &totalTime,
        &mainLoopTime,
        &readControlsTime,
        &setDeltaTTime,
        &pimpleTime,
        &RhoEqnTime,
        &UEqnTime,
        &YEqnTime,
        &EEqnTime,
        &pEqnTime,
        &turbEqnTime,
        &writeTime,
        &rhoFetchTime,
        &YConvectionTime,
        &CombustionModelTime,
        &HeatReleaseTime,
        &YLoopTime,
        &TCEvalTime,
        &ODESolveTime,
        &JacobianEvalTime,
        &dYdTEvalTime,
        &ReactionRateEvalTime,
        &OmegaEvalTime
    });
    //end

    totalTime.start();

    while (runTime.run())
    {
        mainLoopTime.start();

        readControlsTime.start();
        #include "readTimeControls.H"
        readControlsTime.stop();

        setDeltaTTime.start();
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
            // begin copy-pasted YEqn to avoid annoying inclusion of
            // timers
            YConvectionTime.start();
            tmp<fv::convectionScheme<scalar>> mvConvection
            (
                fv::convectionScheme<scalar>::New
                (
                    mesh,
                    fields,
                    phi,
                    mesh.divScheme("div(phi,Yi_h)")
                )
            );
            YConvectionTime.stop();
            {
                CombustionModelTime.start();
                reaction->correct();
                CombustionModelTime.stop();
                HeatReleaseTime.start();
                Qdot = reaction->Qdot();
                HeatReleaseTime.stop();
                SetYInertTime.start();
                volScalarField Yt(0.0*Y[0]);
                SetYInertTime.stop();

                YLoopTime.start();
                forAll(Y, i)
                {
                    if (i != inertIndex && composition.active(i))
                    {
                        volScalarField& Yi = Y[i];

                        fvScalarMatrix YiEqn
                        (
                            fvm::ddt(rho, Yi)
                          + mvConvection->fvmDiv(phi, Yi)
                          - fvm::laplacian(turbulence->muEff(), Yi)
                         ==
                            reaction->R(Yi)
                          + fvOptions(rho, Yi)
                        );

                        YiEqn.relax();

                        fvOptions.constrain(YiEqn);

                        YiEqn.solve(mesh.solver("Yi"));

                        fvOptions.correct(Yi);

                        Yi.max(0.0);
                        Yt += Yi;
                    }
                }
                YLoopTime.stop();

                SetYInertTime.start();
                Y[inertIndex] = scalar(1) - Yt;
                Y[inertIndex].max(0.0);
                SetYInertTime.stop();
            }
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
        rhoFetchTime.stop();

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
