{
    if (npBsubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
/*
        surfaceScalarField alphaPhiSum
        (
            IOobject
            (
                "alphaPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(phi.dimensions(), Zero)
        );
 
*/
        for
        (
            subCycle<volScalarField> pBSubCycle(pB, npBsubCycles);
            !(++pBSubCycle).end();
        )
        {
            #include "pBEqn.H"
//            alphaPhiSum += (runTime.deltaT()/totalDeltaT)*alphaPhi;
        }
 
//        alphaPhi = alphaPhiSum;
    }
    else
    {
        #include "pBEqn.H"
    }
 
}
