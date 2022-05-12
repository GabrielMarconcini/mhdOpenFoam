// --- B-PISO loop
        while (bpiso.correct())
        {
            
            fvVectorMatrix BEqn
            (
                fvm::ddt(B)
              + fvm::div(phiInc, B)
              - fvc::div(phiB, U)
              - fvm::laplacian(DB, B)
            );

            BEqn.solve();

            volScalarField rAB(1.0/BEqn.A());
            surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

            phiB = fvc::flux(B);

            while (bpiso.correctNonOrthogonal())
            {
                fvScalarMatrix pBEqn
                (
                    fvm::laplacian(rABf, pB) == fvc::div(phiB)
                );

                pBEqn.solve(mesh.solver(pB.select(bpiso.finalInnerIter())));

                if (bpiso.finalNonOrthogonalIter())
                {
                    phiB -= pBEqn.flux();
                }
            }

            #include "magneticFieldErr.H"
        }
