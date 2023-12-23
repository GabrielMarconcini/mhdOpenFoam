        fvScalarMatrix pBEqn
        (   
            fvm::ddt(pB)+fvc::div(pB_pos*phiv_pos+pB_neg*phiv_neg) // fvm::div(phiInc,pB)
            +   
            sqr(ch)*fvc::div(bf)
            +
            //sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit)))*pB
            fvm::SuSp(sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit))),pB)
            +0.5*fvm::SuSp(fvc::div(U),pB)
            //+0.5*pB*fvc::div(U) // termo do artigo do Tricco & Price (2012)
            //fvm::Su(sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit)))*pB,pB)
            //sqr(ch/(cd+ROOTVSMALL*sqrt(oneDiffusionUnit)))*pB
        );  

        pBEqn.solve();
        pB_pos = interpolate(pB, pos);
        pB_neg = interpolate(pB, neg);
