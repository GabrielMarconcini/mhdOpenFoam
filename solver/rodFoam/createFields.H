    Info<< "Reading physicalProperties\n" << endl;
    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    dimensionedScalar muMag
    (
    	"muMag",
    	dimensionSet (1, 1, -2, 0, -2, 0, 0),
    	physicalProperties
    );


    Info<< "Reading field sigma\n" << endl;
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField ElPot
    (
        IOobject
        (
    	    "ElPot",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::MUST_READ,
    	    IOobject::AUTO_WRITE
    	),
    	mesh
    );
    
    
    Info<< "Reading field A\n" << endl;
    volVectorField A
    (
        IOobject
        (
    	    "A",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::MUST_READ,
    	    IOobject::AUTO_WRITE
    	),
    	mesh
    );
    
    
    Info<< "Calculating magnetic field B \n" << endl;
    volVectorField B
    (
        IOobject
        (
    	    "B",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
    	),
    	fvc::curl(A)
    );
    
    volVectorField Je
    (
       IOobject
       (
    	    "Je",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::NO_READ,
    	    IOobject::AUTO_WRITE
    	),
    	-sigma*(fvc::grad(ElPot))
    );
