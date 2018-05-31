YConvection.start();
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
YConvection.stop();

{
    CombustionModelTimer.start();
    reaction->correct();
    CombustionModelTimer.stop();
    HeatReleaseTimer.start();
    Qdot = reaction->Qdot();
    HeatReleaseTimer.stop();
    SetYInertTime.start();
    volScalarField Yt(0.0*Y[0]);
    SetYInertTime.stop();

    YLoop.start();
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
