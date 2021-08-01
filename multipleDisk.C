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
    Add volume force in an actuator disk region from thrust, torque 
    and geometry defined in fvSolution

    Use with multipleDiskSimpleFoam

    Modified by Edgar Martínez, August 2020.

    The actuator disk can be defined by adding the following 
    lines in fvSolutions:

    actuatorDiskID
    {
        interiorRadius      1.6;    
        exteriorRadius      20.5;   
        thrust              47.5e3; 
        torque              112.0e3;    
        density             1.2; 
        startPoint          (103.0 0 0);    
        endPoint            (102.0 0 0);   
    }
\*---------------------------------------------------------------------------*/

#include "multipleDisk.H"

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"



namespace Foam {

    defineTypeNameAndDebug(multipleDisk, 0);

        // Default constructor
        multipleDisk::multipleDisk() {
       
        mPointStartCenterLine.x() = 0.0;
        mPointStartCenterLine.y() = 0.0;
        mPointStartCenterLine.z() = 0.0;

        mPointEndCenterLine.x() = 0.0;
        mPointEndCenterLine.y() = 0.0;
        mPointEndCenterLine.z() = 0.0;

        mExtRadius = 0.0;
        mIntRadius = 0.0;
        mThrust = 0.0;
        mTorque = 0.0;
        mRho = 1.0;
        ID = "NULL";
    }



    // Destructor
    multipleDisk::~multipleDisk() {}



    // Convert ID to string
    void multipleDisk::setID(int ID_){
    ID = std::to_string(ID_);
    }



    void multipleDisk::ReadGeometry(const fvMesh &iMesh) {

        if(debug >= 1) { Info << "Reading actuator disk geometry.\n"; }

        std::string object;
        object.append("actuatorDisk");
        object.append(ID);  

        Istream& is1 = iMesh.solutionDict().subDict(object).lookup("interiorRadius");
        is1.format(IOstream::ASCII);
        is1 >> mIntRadius;

        Istream& is2 = iMesh.solutionDict().subDict(object).lookup("exteriorRadius");
        is2.format(IOstream::ASCII);
        is2 >> mExtRadius;

        Istream& is3 = iMesh.solutionDict().subDict(object).lookup("thrust");
        is3.format(IOstream::ASCII);
        is3 >> mThrust;

        Istream& is4 = iMesh.solutionDict().subDict(object).lookup("torque");
        is4.format(IOstream::ASCII);
        is4 >> mTorque;     

        Istream& is5 = iMesh.solutionDict().subDict(object).lookup("density");
        is5.format(IOstream::ASCII);
        is5 >> mRho;

        Istream& is6 = iMesh.solutionDict().subDict(object).lookup("startPoint");
        is6.format(IOstream::ASCII);
        is6 >> mPointStartCenterLine;

        Istream& is7 = iMesh.solutionDict().subDict(object).lookup("endPoint");
        is7.format(IOstream::ASCII);
        is7 >> mPointEndCenterLine;  

        if(debug >= 2) {
            Info << "Actuator disk values loaded from fvSolution:\n";
            Info << "mIntRadius: " << mIntRadius << "\n";
            Info << "mExtRadius: " << mExtRadius << "\n";
            Info << "mThrust: " << mThrust << "\n";
            Info << "mTorque: " << mTorque << "\n";
            Info << "mRho: " << mRho << "\n";
            Info << "mPointStartCenterLine: " << mPointStartCenterLine << "\n";
            Info << "mPointEndCenterLine: " << mPointEndCenterLine << "\n";
        }
    }



    void multipleDisk::CalcActuatorDiskVolForce(const fvMesh &iMesh, vectorField & ioVolumeForce) {

        if(debug >= 1) { Info << "Calculating volume force from actuator disk.\n"; }

        ReadGeometry(iMesh);

        scalar RadialDist2;
        vector LineTangent;
        vector CircumferentialDirection;
        vector TotalForce(0.0, 0.0, 0.0);
        scalar TotalTorque = 0.0;
        scalar DiskVolume = 0;

        // Loop over all cells and check if the cell center is in the actuator disc region
        for (label i = 0; i < iMesh.C().size(); i++) {

            if(PointIsInDisk(mPointStartCenterLine, 
                 mPointEndCenterLine, 
                 iMesh.C()[i],
                         RadialDist2, 
                 LineTangent, 
                 CircumferentialDirection)) {

                if(debug >= 3) { Info << "Point: " << i << " is in the actuator disk. Coordinates: " << iMesh.C()[i] << "\n"; }

                // force is [N / m**3] divided by rho since this is an incompressible solver
                vector axialForce = LineTangent * CalcAxialForce(sqrt(RadialDist2), mRho) / mRho;

                // accumulate force
                ioVolumeForce[i] += axialForce;

                // total axial force, units are [N / rho]
                TotalForce += axialForce*iMesh.V()[i];

                // force is [N / m**3] divided by rho since this is an incompressible solver
                vector circForce = CircumferentialDirection * CalcCircForce(sqrt(RadialDist2), mRho) / mRho;

                // accumulate force
                ioVolumeForce[i] += circForce;

                // total torque, units are [N-m / rho]
                TotalTorque += (CalcCircForce(sqrt(RadialDist2), mRho) / mRho) * sqrt(RadialDist2) * iMesh.V()[i];

                DiskVolume += iMesh.V()[i];
            }
        }

        Info << "Total axial force: " << TotalForce * mRho << " [N]" << "\n";
        Info << "Total torque: " << TotalTorque * mRho << " [N-m]" << "\n";
        Info << "Total disk volume: " << DiskVolume << "\n";
    }



    void multipleDisk::WriteVTK() {

    FILE *file;
    char fileName[100];

    // The cylindrical surface is visualized as 20 rectangular surfaces
    unsigned int NumCells = 20; 

    // 40 points are required ; 20 points at each end of the cylinder
    unsigned int NumPoints = 40; 

    // Number of integers needed in the the VTK file; 
    // each surface has 4 corner points, 
    // so we need 4 corner indices + the index of the surface = 5 indices per surface
    unsigned int NumInts = 5*NumCells; 

    vectorField points(NumPoints, vector::zero);

    // start point is located downwind while end point is located upwind
    vector VecLineTangent(mPointEndCenterLine - mPointStartCenterLine);
    scalar LineTangentLength = sqrt(VecLineTangent.x()*VecLineTangent.x() + VecLineTangent.y()*VecLineTangent.y() + VecLineTangent.z()*VecLineTangent.z());

    if(LineTangentLength!= 0.0) 
    { 
        // normalize the vector
        VecLineTangent /= LineTangentLength; 
    }
    else 
    {
        Info << "Warning: The centerline tangent has zero length.\n";
        return; 
    }

    // We need to find a vector in the radial direction . This can be any vector as long as it points in the radial direction .
    // First try with (1 0 0) and see if we can project it onto the normal plane of the actuator disk resulting in a vector in
    // the radial direction .
    vector VecRadialDirection (1.0, 0.0, 0.0);
    VecRadialDirection -= (VecRadialDirection & VecLineTangent) * VecLineTangent;

    if(mag(VecRadialDirection) < SMALL) {
        // If we enter this if statement , our guess (1 0 0) was parallel to the centerline of the actuator disk . Then
        // we try (0 1 0) instead . Since (1 0 0) was parallel to the centerline, (0 1 0) will for sure not be parallel to
        // the centerline .
        VecRadialDirection.x() = 0.0;
        VecRadialDirection.y() = 1.0;
        VecRadialDirection.z() = 0.0;
        VecRadialDirection -= (VecRadialDirection & VecLineTangent) * VecLineTangent;
    }

    if(mag(VecRadialDirection) > SMALL) {
        VecRadialDirection /= mag(VecRadialDirection);
    }
    else {
        Info << "Warning in actuatorDiskExplicitForce::WriteVTK(): mag(VecRadialDirection) close to zero.\n";
    }

    vector VecRadialDirection2 = VecLineTangent ^ VecRadialDirection;
    scalar XLocal = 0.0, YLocal = 0.0;

    // Compute points on first side of disk region
    double phi = 0.0;
    for(unsigned int i = 0; i < NumCells; i ++) {
        XLocal = mExtRadius*cos(phi);
        YLocal = mExtRadius*sin(phi);
        vector point(mPointStartCenterLine + XLocal*VecRadialDirection + YLocal*VecRadialDirection2);
        points[i] = point;
        phi += (1.0/ double(NumCells)) * 2 *mPI;
    }

    // Compute points on second side of disk region
    phi = 0.0;
    for(unsigned int i = 0; i < NumCells; i++) {
        XLocal = mExtRadius*cos(phi);
        YLocal = mExtRadius*sin(phi);
        vector point(mPointEndCenterLine + XLocal*VecRadialDirection + YLocal*VecRadialDirection2);
        points[NumCells + i] = point;
        phi += (1.0/ double(NumCells))*2*mPI;
    }

    std::string s;
    s.append("actuatorDisk");
    s.append(ID);
    s.append(".vtk");
    sprintf(fileName, s.c_str());
    file = fopen(fileName,"w");

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Analytical surface of actuator disk.\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(file, "POINTS %i float\n", NumPoints);

    for(int i = 0; i < points.size(); i++) {
        fprintf(file, "%e %e %e\n", points[i].x(), points[i].y(), points[i].z()) ;
    }

    fprintf(file, "CELLS %i %i\n", NumCells, NumInts);
    
    for(unsigned int i = 0; i < NumCells-1; i++) {
        fprintf(file, "%i %i %i %i %i\n",4,i,i+NumCells,i+NumCells+1,i+1);
    }
    
    fprintf(file,"%i %i %i %i %i\n",4,NumCells-1,2*NumCells-1,NumCells,0);

    fprintf(file,"CELL_TYPES %i\n", NumCells);

    for(unsigned int i = 0; i < NumCells; i++) {
        fprintf(file, "%i\n",9);
    }
    
    fclose(file);
    }



    bool multipleDisk::PointIsInDisk(const vector &iPointStartCenterLine, 
                    const vector &iPointEndCenterLine,
                    const vector &iPoint, 
                    scalar &oDist2, 
                    vector &oLineTangent, 
                    vector &oCircumferentialDirection) {

        // center line vector "v":
        vector VecLineTangent(iPointEndCenterLine - iPointStartCenterLine);

        // cente line vector magnitude:
        scalar LineTangentLength = sqrt(VecLineTangent.x()*VecLineTangent.x() + 
                            VecLineTangent.y()*VecLineTangent.y() +
                            VecLineTangent.z()*VecLineTangent.z());

        // normalize if length != 0:
        if(LineTangentLength != 0.0) {
            VecLineTangent /= LineTangentLength;
        }
        else
        {
            Info << "Warning: The centerline tangent has zero length.\n";
            return false;
        }

        // pass center line vector to argument:
        oLineTangent = VecLineTangent;

        // start point to test point vector "s = q - ps"
        vector VecStartLineToPoint(iPoint - iPointStartCenterLine);

        // projection of "s" onto "v"
        scalar PointProjOnLine = VecStartLineToPoint & VecLineTangent;

        // check if the point is inside the actuator disk in the axial direction
        // is d > 0 and d > L?
        if(!(PointProjOnLine >= 0.0 && PointProjOnLine <= LineTangentLength)) {
            return false;
        }
    
        // radial vector "r"
        vector VecLineToPoint(VecStartLineToPoint - (VecLineTangent * PointProjOnLine));

        // squared magnitude of "r"
        scalar RadialDist2 = VecLineToPoint.x()*VecLineToPoint.x() + 
                             VecLineToPoint.y()*VecLineToPoint.y() +
                             VecLineToPoint.z()*VecLineToPoint.z();

        // pass squared magnitude of "r" to argument:
        oDist2 = RadialDist2;

        // tangent vector:
        oCircumferentialDirection = VecLineTangent ^ VecLineToPoint;

        // normalize tangent vector:
        oCircumferentialDirection /= mag(oCircumferentialDirection);

        // check mag(r) < Rp && mag(r) > Rh
        return (RadialDist2 <= mExtRadius*mExtRadius && RadialDist2 >= mIntRadius*mIntRadius);
    }



    scalar multipleDisk::CalcAxialForce(const scalar &iRadialDist, const scalar &iRho) {

    // alias
    scalar pi = mPI;
    scalar r = iRadialDist/mExtRadius;
    scalar rh = mIntRadius/mExtRadius;
    scalar ri = mIntRadius;
    scalar ro = mExtRadius;
    scalar T = mThrust;
    scalar D = CalcDiskThickness();

    // auxiliary
        scalar axialForce = 0.0;
        scalar rast = (r - rh) / (1.0 - rh);

    // computations
        scalar Ax = (105.0 / 8.0) * T / (D * pi * (3.0 * ri + 4.0 * ro) * (ro - ri));
        axialForce = Ax * rast * sqrt(1.0 - rast);

        if (debug >= 2) {
            Info << "Axial force: " << axialForce << "\n";
        }

        return axialForce;
    }



    scalar multipleDisk::CalcCircForce(const scalar &iRadialDist, const scalar &iRho) {

    // alias
    scalar pi = mPI;
    scalar r = iRadialDist/mExtRadius;
    scalar rh = mIntRadius/mExtRadius;
    scalar ri = mIntRadius;
    scalar ro = mExtRadius;
    scalar Q = mTorque;
    scalar D = CalcDiskThickness();

    // auxiliary
        scalar tangentialForce = 0.0;
        scalar rast = (r - rh) / (1.0 - rh);

    // computations
        scalar At = (105.0 / 8.0) * Q / (D * pi * ro * (ro - ri) * (3.0 * ro + 4.0 * ri));
        tangentialForce = (At* rast * sqrt(1.0 - rast) / (rast * (1.0 - rh) + rh));

        if(debug >= 2) {
            Info << "Tangential force: " << tangentialForce << "\n";
        }

        return tangentialForce;
    }
} 
