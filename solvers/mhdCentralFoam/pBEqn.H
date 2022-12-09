        fvScalarMatrix pBEqn
        (   
            fvm::ddt(pB) 
            +   
            sqr(ch)*fvc::div(bf)
            +   
            fvm::SuSp(sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit))),pB)
            //fvm::Su(sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit)))*pB,pB)
            //sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit)))*pB
        );  

        pBEqn.solve();
        pB.correctBoundaryConditions();
    
        #include "magneticFieldErr.H" 
