---
title: "Creating your own OpenFOAM solver"
date: 2022-02-08T19:51:03+05:30

categories: [OpenFOAM]
tags: [OpenFOAM development]
toc: false
author: "Satya Sasank"
---
AIM:

The aim of this project is to create a custom solver in OpenFOAM that solves a transport equation for a scalar named s.

OBJECTIVES:

The solver must be modified from icoFoam which is present in $FOAM/applications/solvers/incompressible/icoFoam.
The case should be a scalar quantity 's' that is being transported by velocity u and the transport equation is represented by `{del u }/{del t } + nabla. (us) = 0`
The case must be modified from $FOAM_TUTORIALS/incompressible/icoFoam/cavity.

PROCEDURE:

$FOAM/applications/solvers/incompressible/icoFoam and $FOAM_TUTORIALS/incompressible/icoFoam/cavity to this directory. We rename these folders to scalarFoam and scalarCavity. Our objective here is to add the scalar transport equation to the solver. So we first start by making changes in the solver.

scalarFoam.c

we first go into the scalarFoam and rename the icoFoam.c to scalarFoam.c. Now in the scalarFoam.c we need to add the scalar transport equation of s. we can do this by adding the following piece of code at the end of the PISO loop.

```c++
fvScalarMatrix sEqn
        {
            fvm::ddt(s)
            + fvm::div(phi,s)
            
        };
        sEqn.relax();
        sEqn.solve();
```

In the above piece of code, we are solving for a scalar 's' and the solution is to be stored in a fvScalarMatrix called sEqn. Inside this, we are defining the transport equation. The first term indicates the time derivative and the second term indicates the convection term.

Finally, we add the necessary relaxations and solve the equation. In the above equation, we have defined a new variable called s. We need to define this variable, We can do this by creating a file for s in the 0 folder of the case directory.

S
```c++
/*--------------------------------*- C++ -*----------------------------------*
| =========                 |                                                 |
|       /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|      /   O peration     | Version:  v2012                                 |
|     /    A nd           | Website:  www.openfoam.com                      |
|    /     M anipulation  |                                                 |
*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];


internalField   uniform 0;

boundaryField
{
    movingWall
    {
        type            zeroGradient;
    }
    fixedWalls
    {
        type            zeroGradient;
    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
```




The above file is created by copying the pressure file and making the necessary changes. Here the dimensions of the scalar can be anything but for this project, the same dimensions of time are used. Now we need to tell OpenFOAM that s is defined here. For that, we add the following piece of code to the createFields.H

createFields.H
```c++
Info<< "Reading field sn" << endl;
volScalarField s
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
```
Now our solver is ready we need to compile it for that we make the following changes.

We go into the make folder and files and modify them as shown

files

scalarFoam.C
```c++
EXE = $(FOAM_USER_APPBIN)/scalarFoam
```

The options file must not be changed since it contains the libraries required. Once this is done we open a terminal in the scalarFoam folder and run wmake. This will compile the solver.

Now we must tell OpenFoam what kind of schemes to use for solving s. For that we goto fvSchemes

fvSchemes
```c++
divSchemes
{
    default         none;
    div(phi,U)      Gauss linear;
    div(phi,s)      Gauss upwind;
}
```
Now we need to specify what kind of solver to be used for solving s. 

fvSolution
```c++
s
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
```

We cant use the same exact solver we used for solving p. This is because PCG is meant for symmetrical, PBiCG is meant for asymmetrical matrices, and smoothsolver for either of those.

We shall now specify a non-uniform initial condition for the s. The default field values set the default value of the fields, i.e. the value the field takes unless specified otherwise in the regions sub-dictionary. That sub-dictionary contains a list of subdictionaries containing field values that override the defaults in a specified region. The region creates a set of points, cells, or faces based on some topological constraint. Here, boxToCell creates a bounding box within a vector minimum and maximum to define the set of cells of the water region. This will be done by running the set fields utility. It requires a setFieldsDict dictionary, located in the system directory, whose entries for this case are shown below.

The following piece of code will be stored in a file called setFields and will be stored in the system folder.
```c++
/*--------------------------------*- C++ -*----------------------------------*
| =========                 |                                                 |
|       /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|      /   O peration     | Version:  v2012                                 |
|     /    A nd           | Website:  www.openfoam.com                      |
|    /     M anipulation  |                                                 |
*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      setFieldsDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defaultFieldValues
(
    volScalarFieldValue s 0
);

regions
(
    sphereToCell
    {
        centre (0.05 0.05 0);
        radius 0.03; 

        fieldValues
        (
            volScalarFieldValue s 1
        );
    }
);


// ************************************************************************* //
```

This is taken from the dambreak problem and is configured for the scalar s

Once this is done we are ready to run the simulation

Before that, we need to make some changes in the time step and mesh size 

In the blockmeshdict file, we set the grading as follows
```c++
blocks
(
    hex (0 1 2 3 4 5 6 7) (100 100 1) simpleGrading (1 1 1)
);
```

 Also in the controlDict file, the time size is given as follows. 
```c++
application     scalarFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         2;

deltaT          0.0005;
```
The objective of this simulation is to observe the change in the scalar field s with time, in the cavity setup.

Now that everything is setup we execute the following commands

blockMesh

checkMesh

setFields

scalarFoam

paraFoam

This will open ParaView and the results will be as shown.

RESULTS:

at t = 0
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634920915.jpg)

at t = 0.1s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921112.jpg)

at t = 0.2s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921128.jpg)

at t = 0.3s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921143.jpg)

at t = 0.4s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921167.jpg)

at t = 0.5s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921182.jpg)

at t = 0.6s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921204.jpg)

at t = 0.7s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921215.jpg)

at t = 0.8s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921234.jpg)

at t = 0.9s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921250.jpg)

at t = 1s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921268.jpg)

at t = 2s
![](https://sklc-tinymce-2021.s3.amazonaws.com/comp/2021/10/ScalarCavity_1634921281.jpg)

CONCLUSION:

In the final time step, we can see that how the scalar quantity has dispersed in the domain due to the presence of velocity. Since we have considered only the convection term and not the diffusion term we cannot observe the diffusion of the scalar in our domain.


