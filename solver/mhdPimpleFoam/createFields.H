#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

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

#include "compressibleCreatePhi.H"
#include "createIncompressiblePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

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

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"


//--- Variáveis de campos magnéticas:
Info<< "Reading transportProperties\n" << endl;

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

//dimensionedScalar nu
//(
//    "nu",
//    dimViscosity,
//    transportProperties
//);

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

//--- Lendo campo de pressão magnética:
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

//--- Criando campo de fluxo magnético:
#include "createPhiB.H"

dimensionedScalar DB = 1.0/(mu0*sigma);
DB.name() = "DB";

dimensionedScalar DBU = 1.0/(2.0*mu0);
DBU.name() = "DBU";

//--- Constantes magnéticas da equação da energia:
dimensionedScalar EDB = 1.0/(mu0*mu0*sigma);
EDB.name() = "EDB";

dimensionedScalar EDBU = 1.0/(2.0*mu0*mu0*sigma);
EDBU.name() = "EDBU";



//label pRefCell = 0;
//scalar pRefValue = 0.0;
//setRefCell(p, piso.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(pB.name());
