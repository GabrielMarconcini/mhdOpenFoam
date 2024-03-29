{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K)  + fvc::div(phi, K)
      + fvc::ddt(KB)      + fvc::div(phiInc,KB)
      + fvc::div(phiInc,KB)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
      - fvc::div(phiB, 2.0*DBU*(B & U))
      - fvc::laplacian(EDBU*magSqr(B))
      + fvc::div
        (
            fvc::div(phiB, EDB*B), 
            "div(div(BB))"
        )
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
