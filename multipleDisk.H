/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    multipleDisk

Description
    Add volume force in an actuator disk region from thrust, torque and 
    geometry defined in fvSolution

    Use with multipleDiskSimpleFoam

    Modified by Edgar Martínez, august 2020.
\*---------------------------------------------------------------------------*/

#ifndef MULTIPLEDISK_H_
#define MULTIPLEDISK_H_

#include "fvCFD.H"

namespace Foam {

    class multipleDisk {

    public:

    // Runtime type information 
    TypeName("multipleDisk");

    // Constructor
    multipleDisk();

    // Destructor
    ~multipleDisk();

    // Reads dictionary input from the fvSolution file
    void ReadGeometry(const fvMesh &iMesh);

    // Axial and tangential forces applied on each cell of the disc
    void CalcActuatorDiskVolForce(const fvMesh &iMesh, vectorField &ioVolumeForce);

    // Write the VTK file, this can be visualized in Paraview
    void WriteVTK();

    // Set the turbine ID, or number of turbine
    void setID(int ID_);

    private:

    // Start point vector
    vector mPointStartCenterLine;

    // End point vector
    vector mPointEndCenterLine;

    // External radius
    scalar mExtRadius;

    // Internal radius
    scalar mIntRadius;

    // Prescribed thrust [N]
    scalar mThrust;

    // Prescribed torque [N-m]
    scalar mTorque;

    // Fluid density
    scalar mRho;

    // ID string. ID argument is an int.
    std::string ID;

    // Pi constant
    static constexpr double mPI = 3.141592654;

    // Check whether the point is inside the disc
    bool PointIsInDisk(const vector &iPointStartCenterLine, 
            const vector &iPointEndCenterLine, 
            const vector &iPoint, 
            scalar &oDist2, 
            vector &oLineTangent, 
            vector &oCircumferentialDirection);
    
    // Calculate axial force scalar 
    scalar CalcAxialForce(const scalar &iRadialDist, const scalar &iRho);

    // Calculate tangential force scalar
    scalar CalcCircForce(const scalar &iRadialDist, const scalar &iRho);

    // Calculate disc thickness
    scalar CalcDiskThickness() { return mag(mPointEndCenterLine - mPointStartCenterLine); };
  };
}
#endif
