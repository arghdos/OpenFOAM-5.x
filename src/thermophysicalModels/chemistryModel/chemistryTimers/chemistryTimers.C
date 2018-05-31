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

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "StopWatch.H"
#include "CustomTimers.H"

// create timers
namespace Foam
{
    StopWatch QDotEvalTime(string("heat-release rate evaluation"));
    //- chemical time scale eval
    StopWatch TCEvalTime(string("chemical time-scale evaluation"));
    // ODE-solver timers
    StopWatch ODESolveTime(string("ode integration"));
    StopWatch JacobianEvalTime(string("jacobian evaluation"));
    StopWatch dYdTEvalTime(string("dydt evaluation"));
    StopWatch ReactionRateEvalTime(string("reaction rate evaluation"));
    StopWatch OmegaEvalTime(string("species rate evaluation"));
}
