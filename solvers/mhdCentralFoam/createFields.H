#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& e = thermo.he();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

//==========================================================================================================================================
//--- Variáveis magnéticas (adicionado):

Info<< "Reading transportProperties\n" <<endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar mu0
(
    "mu0",
    dimensionSet(1, 1, -2, 0, 0, -2, 0),
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
);

dimensionedScalar Cr
(
    "Cr",
    dimless, 
    transportProperties
);

//--- Lendo pB:
Info<< "Reading field pB\n" << endl;

volScalarField pB
(
    IOobject
    (
        "pB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//--- Lendo campo magnético:
Info<< "Reading field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField  delta
(
   IOobject
   (
       "delta",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("delta", dimless, 0.0)
);

scalar deltaMax=ROOTVSMALL;
const surfaceScalarField& deltaCoeffs=mesh.deltaCoeffs();


forAll(mesh.C(), celli)
{
      deltaMax = ROOTVSMALL;
      const labelList& faces = mesh.cells()[celli];
      forAll(faces,facei)
      {
            deltaMax=max(deltaCoeffs[facei],deltaMax);
      }
      delta.primitiveFieldRef()[celli]=deltaMax;
}

volScalarField ch
(
     IOobject
     (   
         "ch",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),  
     mesh,
     dimensionedScalar("ch", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
);  


volScalarField cd
(   
     IOobject
     (   
         "cd",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("cd", dimensionSet(0, 1, -0.5, 0, 0, 0, 0), 0.0)
);


//--- Constantes magnéticas das equações:
dimensionedScalar DB = 1.0/(mu0*sigma);
DB.name() = "DB";

dimensionedScalar DBU = 1.0/(2.0*mu0);
DBU.name() = "DBU";

dimensionedScalar EDB = 1.0/(mu0*mu0*sigma);
EDB.name() = "EDB";

dimensionedScalar EDBU = 1.0/(2.0*mu0*mu0*sigma);
EDBU.name() = "EDBU";

mesh.setFluxRequired(pB.name());
//============================================================================================================================================

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U)) + DBU*magSqr(B)
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

//#include "createIncompressiblePhi.H"  //adicionado

