/*
 *  Copyright (c) 2011       Marius Cautun
 *
 *                           Kapteyn Astronomical Institute
 *                           University of Groningen, the Netherlands
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*!
Here are the functions used for data input and output.
To easily find how to read the input data go to the function:
    


*/

#ifndef INPUT_OUTPUT_HEADER
#define INPUT_OUTPUT_HEADER

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
namespace bfs=boost::filesystem;

#include "define.h"
#include "particle_data.h"
#include "user_options.h"
#include "message.h"


// contains the definitions of some classes and functions used only for input and output purposes
#include "io/input_output.h"

// different data format readers and writters
#include "io/gadget_reader.cc"  // reader for Gadget snapshots
#include "io/binary_io.cc"      // reader/writer for binary file format
#include "io/text_io.cc"        // reader/writer for text files
#include "io/my_io.cc"          // reader/writer for your custom files






//! Functions for reading the input data

/* In the next few lines you can easily select which function to use to read the input data. This is done via the precompiler variable 'DATA_INPUT_FUNCTION'. */

/* Read the data from a text file - see the "src/io/text_io.cc" for the definition of this function */
#define DATA_INPUT_FUNCTION readTextFile

/* Read only positions from a text file - see the "src/io/text_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readTextFile_positions

/* Read only positions from a text file and on top of that also reads user defined sampling points from an additional file (USE ONLY WHEN YOU NEED CUSTOM SAMPLING POINTS) - see the "src/io/text_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readTextFile_userDefinedSampling

/* Read the input data from a Gadget snapshot file - see the "src/io/gadget_reader.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readGadgetFile

//#define DATA_INPUT_FUNCTION readMyFile

/* Read the input data from a binary file - see the "src/io/binary_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readBinaryFile

/* Read the input data from your custom type file - see the "src/io/my_io.cc" for the definition of this function */
//#define DATA_INPUT_FUNCTION readMyFile

/* This function reads the particle data used for the DTFE computation. You can modify this function as you please. */
void readInputData(std::vector<Particle_data> *p,
                   std::vector<Sample_point> *samplingCoordinates,
                   User_options *userOptions)
{
    std::string filename = userOptions->inputFilename;  // the name of the input file
    Read_data<Real> readData; // will store the data read from the input file
    
    
    // Read the data from the input file - see the definition of 'DATA_INPUT_FUNCTION' before the beginning of this function to select a different type of input data file
    DATA_INPUT_FUNCTION( filename, &readData, userOptions );
    
    
    // 'userOptions->MpcValue' is the conversion factor from the units in the input data file to Mpc units - do the next computation only if userOptions->MpcValue!=1
    if ( userOptions->MpcValue!=Real(1.) )
    {
        for (size_t i=0; i<userOptions->boxCoordinates.size(); ++i)
            userOptions->boxCoordinates[i] /= userOptions->MpcValue;
        
        size_t noParticles = readData.noParticles();  // returns the number of particles
        float *positions = readData.position();  //returns pointer to array storing particle positions
        for (size_t i=0; i<noParticles*NO_DIM; ++i)
            positions[i] /= userOptions->MpcValue;
        
        // if the user inserted user defined sampling points, divide the coordinates of those points by the normalization factor
        if ( readData.noSamples()!=size_t(0) )
        {
            float *sampling = readData.sampling();  // returns pointer to user defined sampling coordinates
            float *delta = readData.delta();        // returns pointer to user defined sampling cell sizes
            for (size_t i=0; i<readData.noSamples()*NO_DIM; ++i)
            {
                sampling[i] /= userOptions->MpcValue;
                delta[i] /= userOptions->MpcValue;
            }
        }
    }
    
    
    // now store the data in the 'Particle_data list'. It also copies the user given sampling coordinates, if any - none in this case.
    readData.transferData( p, samplingCoordinates );
}

//! Functions for writing the output data

/* In the next few lines you can easily select which function to use to write the output data to file. This is done via the precompiler variable 'DATA_OUTPUT_FUNCTION'. */
/* Writes the data to a binary file - see "binary_io.cc" for function definition */
//#define DATA_OUTPUT_FUNCTION  writeBinaryFile

//#define DATA_OUTPUT_FUNCTION  writeMyFile

/* Writes the data to a text file - see "text_io.cc" for function definition */
//#define DATA_OUTPUT_FUNCTION  writeTextFile

/* Writes the data to a text file, but on each line writes also the indices of the grid cell corresponding to the result being written - see "text_io.cc" for function definition */
//#define DATA_OUTPUT_FUNCTION  writeTextFile_gridIndex

/* Writes the data to a text file, but on each line writes also the sampling point coordinates corresponding to the result being written - see "text_io.cc" for function definition */
#define DATA_OUTPUT_FUNCTION  writeTextFile_samplingPosition

/* Writes the data to a text file, but on each line writes also the sampling point coordinates for a redshift cone grid - see "text_io.cc" for function definition */
//#define DATA_OUTPUT_FUNCTION  writeTextFile_redshiftConePosition

/* Writes the data to your custom type file - write your own custom output function in file "my_io.cc" */
//#define DATA_OUTPUT_FUNCTION  writeMyFile

/* This function writes the output data to a file, each different 'variable' being written to a separate file. 
You can modify this function as you please. */
void writeOutputData(Quantities &uQuantities,
                     Quantities &aQuantities,
                     User_options const &userOptions)
{
    // output the desired quantities to file/files
    // outputs the density
    if ( userOptions.uField.density )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.density, userOptions.outputFilename + ".den", "density", userOptions );
    }
    if ( userOptions.aField.density )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.density, userOptions.outputFilename + ".a_den", "volume averaged density", userOptions );
    }
    
    
    // outputs the velocity
    if ( userOptions.uField.velocity )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.velocity, userOptions.outputFilename + ".vel", "velocity", userOptions );
    }
    if ( userOptions.aField.velocity )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.velocity, userOptions.outputFilename + ".a_vel", "volume averaged velocity", userOptions );
    }
    
    
    // outputs the velocity gradient
    if ( userOptions.uField.velocity_gradient )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.velocity_gradient, userOptions.outputFilename + ".velGrad", "velocity gradient", userOptions );
    }
    if ( userOptions.aField.velocity_gradient )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.velocity_gradient, userOptions.outputFilename + ".a_velGrad", "volume averaged velocity gradient", userOptions );
    }
    
    
    // outputs the velocity divergence
    if ( userOptions.uField.velocity_divergence )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.velocity_divergence, userOptions.outputFilename + ".velDiv", "velocity divergence", userOptions );
    }
    if ( userOptions.aField.velocity_divergence )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.velocity_divergence, userOptions.outputFilename + ".a_velDiv", "volume averaged velocity divergence", userOptions );
    }
    
    
    // outputs the velocity shear
    if ( userOptions.uField.velocity_shear )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.velocity_shear, userOptions.outputFilename + ".velShear", "velocity shear", userOptions );
    }
    if ( userOptions.aField.velocity_shear )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.velocity_shear, userOptions.outputFilename + ".a_velShear", "volume averaged velocity shear", userOptions );
    }
    
    
    // outputs the velocity vorticity
    if ( userOptions.uField.velocity_vorticity )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.velocity_vorticity, userOptions.outputFilename + ".velVort", "velocity vorticity", userOptions );
    }
    if ( userOptions.aField.velocity_vorticity )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.velocity_vorticity, userOptions.outputFilename + ".a_velVort", "volume averaged velocity vorticity", userOptions );
    }
    
    
    // outputs the scalar fields
    if ( userOptions.uField.scalar )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.scalar, userOptions.outputFilename + ".scalar", "scalar", userOptions );
    }
    if ( userOptions.aField.scalar )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.scalar, userOptions.outputFilename + ".a_scalar", "volume averaged scalar", userOptions );
    }
    
    
    // outputs the scalar fields gradient
    if ( userOptions.uField.scalar_gradient )
    {
        DATA_OUTPUT_FUNCTION( uQuantities.scalar_gradient, userOptions.outputFilename + ".scalarGrad", "scalar gradient", userOptions );
    }
    if ( userOptions.aField.scalar_gradient )
    {
        DATA_OUTPUT_FUNCTION( aQuantities.scalar_gradient, userOptions.outputFilename + ".a_scalarGrad", "volume averaged scalar gradient", userOptions );
    }
}























#endif
