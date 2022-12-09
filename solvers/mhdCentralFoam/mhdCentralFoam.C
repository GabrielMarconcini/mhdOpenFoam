/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    mhdCentralFoam

Group
    mhdCompressibleSolvers

Description
    Density-based compressible flow solver based on central-upwind
    schemes of Kurganov and Tadmor for MHD flows.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based solver based on central-upwind schemes of"
        " Kurganov and Tadmor, for MHD compressible flows."
    );

    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "createFvOptions.H"
    

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);
    dimensionedScalar chGlobal(dimLength/dimTime,Zero);    
    dimensionedScalar cdGlobal(dimLength/sqrt(dimTime),Zero);    

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField B_pos(interpolate(B, pos));                    //adicionado: B+
        surfaceVectorField B_neg(interpolate(B, neg));                    //adicionado: B-
        
        surfaceScalarField Bf_pos("Bf_pos", B_pos & mesh.Sf());           //adicionado: Bf+
        surfaceScalarField Bf_neg("Bf_neg", B_neg & mesh.Sf());           //adicionado: Bf-
        

        surfaceScalarField Bn_pos("Bn_pos", Bf_pos/mesh.magSf());         //adicioando: magnetic field component normal to cell face +
        surfaceScalarField Bn_neg("Bn_neg", Bf_neg/mesh.magSf());         //adicionado: magnetic field component normal to cell face -

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);
        
        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);
        
        surfaceScalarField pM_pos("pM_pos", DBU*magSqr(B_pos));           //adicionado: pressão magnética +
        surfaceScalarField pM_neg("pM_neg", DBU*magSqr(B_neg));           //adicionado: pressão magnética -
        
        surfaceScalarField pG_pos("pG_pos", p_pos + pM_pos);              //adicionado: pressão total +
        surfaceScalarField pG_neg("pG_neg", p_neg + pM_neg);              //adicionado: pressão total -
        
        surfaceScalarField pB_pos(interpolate(pB, pos));                  //adicionado: variável potêncial +
        surfaceScalarField pB_neg(interpolate(pB, neg));                  //adicionado: variável potêncial -

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        // Note: extracted out the orientation so becomes unoriented
        phiv_pos.setOriented(false);
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
        phiv_neg.setOriented(false);

        surfaceScalarField bf("bf", 0.5*(Bf_neg + Bf_pos));                   //adicionado: coeficiente bf
 
        volScalarField vS("vS", sqrt(thermo.Cp()/thermo.Cv()*rPsi));          //adicionado: velocidade do som
        surfaceScalarField vS_pos("vS_pos", interpolate(vS, pos, T.name()));  //adicionado: velocidade do som +
        surfaceScalarField vS_neg("vS_neg", interpolate(vS, neg, T.name()));  //adicionado: velocidade do som -
        
        surfaceScalarField vA_pos("vA_pos", mag(B_pos)/(sqrt(mu0*rho_pos)));  //adicionado: velocidade de Alfvén +
        surfaceScalarField vA_neg("vA_neg", mag(B_neg)/(sqrt(mu0*rho_neg)));  //adicionado: velocidade de Alfvén -
        
        surfaceScalarField c_pos
        (
            "c_pos",
            sqrt(0.5*(sqr(vS_pos) + sqr(vA_pos) + sqrt(sqr(sqr(vS_pos) + sqr(vA_pos)) - 4*sqr(vS_pos)*(magSqr(Bn_pos)/(mu0*rho_pos)))))
        );//adicionado: c+ (Equação 24)
        
        surfaceScalarField c_neg
        (
            "c_neg",
            sqrt(0.5*(sqr(vS_neg) + sqr(vA_neg) + sqrt(sqr(sqr(vS_neg) + sqr(vA_neg)) - 4*sqr(vS_neg)*(magSqr(Bn_neg)/(mu0*rho_neg)))))
        );//adicionado: c- (Equação 24)
        
        surfaceScalarField cf("cf", min(c_pos,c_neg));                         //adicionado: cf (Equação 23)

        surfaceScalarField cSf("cSf", cf*mesh.magSf());                        //adicioando: cf*|Sf| (Equação 21 e 22)
        
        
        //volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        //surfaceScalarField cSf_pos
        //(
        //    "cSf_pos",
        //    interpolate(c, pos, T.name())*mesh.magSf()
        //);

        //surfaceScalarField cSf_neg
        //(
        //    "cSf_neg",
        //    interpolate(c, neg, T.name())*mesh.magSf()
        //);

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf, phiv_neg + cSf), v_zero)
            //max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        ); //modificado: Equação 21

        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf, phiv_neg - cSf), v_zero)
            //min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        ); //modificado: Equação 22

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Calculate ch, cd fields and print global values
        scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());
        ch.primitiveFieldRef() = 0.5*(sumAmaxSf/mesh.V().field())/delta;
        cd.primitiveFieldRef() = sqrt((-runTime.deltaTValue()*sqr(ch))/log(Cr));
  
        chGlobal=min(ch.v());
        cdGlobal=sqrt((-runTime.time().deltaT()*sqr(chGlobal))/log(Cr));

        Info<< "chGlobal = " << chGlobal << nl 
            << "cdGlobal = " << cdGlobal << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;
        
        surfaceVectorField phiB(aphiv_pos*B_pos + aphiv_neg*B_neg);                //adicionado: (a*phif+*B+) + ((1-a)*phif-*B-) + wf*(B- + B+)
        phiB.setOriented(true);

        surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);
        // Note: reassembled orientation from the pos and neg parts so becomes
        // oriented
        phiU.setOriented(true);

        surfaceVectorField phiUp(phiU + (a_pos*pG_pos + a_neg*pG_neg)*mesh.Sf() - DBU*((B_neg*bf) + (B_pos*bf))); //modificado

        surfaceScalarField UB_pos("UB_pos", U_pos & B_pos);
        surfaceScalarField UB_neg("UB_neg", U_neg & B_neg);

        UB_pos.setOriented(false);
        UB_neg.setOriented(false);

        surfaceScalarField gradMagSqrB("gradMagSqrB", fvc::snGrad(magSqr(B)));
        gradMagSqrB.setOriented(false);

        surfaceScalarField BdotGradB_pos("BdotGradB_pos", B_pos & fvc::snGrad(B));
        surfaceScalarField BdotGradB_neg("BdotGradB_neg", B_neg & fvc::snGrad(B));
        BdotGradB_pos.setOriented(false);
        BdotGradB_neg.setOriented(false);

        bf.setOriented(false);        

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + DBU*magSqr(B_pos) + pG_pos) - 2.0*a_pos*DBU*bf*UB_pos    //modificado
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + DBU*magSqr(B_neg) + pG_neg) - 2.0*a_neg*DBU*bf*UB_neg    //modificado
          + aSf*pG_pos - aSf*pG_neg
        );

        if (!idealMHD)
        {
           phiEp -= EDBU*gradMagSqrB * mesh.magSf() - EDB*a_pos*BdotGradB_pos * mesh.magSf() - EDB*a_neg*BdotGradB_neg * mesh.magSf();// resistive MHD
        }
       
        bf.setOriented(true);        
 
        surfaceVectorField phiBU
        (
            "phiBU",
            phiB - ((a_pos*bf*U_pos) + (a_neg*bf*U_neg)) + 0.5*(pB_neg + pB_pos)*mesh.Sf()
        );  //adicionado: div(UB - BU) + grad(pB)
        
        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));
        
        // --- Solve induction
        if (idealMHD)
        {
           solve(fvm::ddt(B) + fvc::div(phiBU));                          // ideal MHD
        }
        else
        {
           solve(fvm::ddt(B) + fvc::div(phiBU) - fvm::laplacian(DB, B)); // resistive MHD
        }

 
        B.correctBoundaryConditions();

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() =
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );

       fvScalarMatrix EEqnRCF
        (   
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
            ==  
            fvOptions(rhoE)
        );  

        fvOptions.constrain(EEqnRCF);
        EEqnRCF.solve();

//        solve
//        (
//            fvm::ddt(rhoE) + fvc::div(phiEp)
//          - fvc::div(sigmaDotU)
//        );

        e = rhoE/rho - 0.5*magSqr(U) - (DBU/rho)*magSqr(B);            //modificado
        fvOptions.correct(e);         //Added
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            ) + (DBU.value()*magSqr(B.boundaryField()));               //modificado

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U)) + DBU*magSqr(B);            //modificado
        }

        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
        
        const dimensionedScalar oneDiffusionUnit(dimensionSet(0, 2, -1, 0, 0, 0, 0 ), 1.0 );        

        #include "pBEqn.H"
        
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
