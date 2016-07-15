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



/* This file defines a reader for binary Gadget sanpshot files. The reader given here reads ONLY dark matter particle positions and properties. */



// Header structure for reading Gadget snapshots
struct Gadget_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
	
	// Function that prints the Gadget header.
	void print()
	{
		std::cout << "\nThe header of the Gadget file contains the following info:\n" 
				<< "npart[6]     =  " << npart[0] << "  " << npart[1] << "  " << npart[2] << "  " << npart[3] << "  " <<  npart[4] << "  " <<  npart[5] << "\n"
				<< "mass[6]      =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
				<< "time         =  " << time << "\n"
				<< "redshift     =  " << redshift << "\n"
				<< "flag_sfr     =  " << flag_sfr << "\n"
				<< "flag_feedback=  " << flag_feedback << "\n"
				<< "npartTotal[6]=  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "  " << "\n"
				<< "flag_cooling =  " << flag_cooling << "\n"
				<< "num_files    =  " << num_files << "\n"
				<< "BoxSize      =  " << BoxSize << "\n"
				<< "Omega0       =  " << Omega0 << "\n"
				<< "OmegaLambda  =  " << OmegaLambda << "\n"
				<< "h            =  " << HubbleParam << "\n\n";
	}
	
	// Swap endianness
	void swapBytes()
	{
		ByteSwapArray( npart, 6 );
		ByteSwapArray( mass, 6 );
		BYTESWAP( time );
		BYTESWAP( redshift );
		BYTESWAP( flag_sfr );
		BYTESWAP( flag_feedback );
		ByteSwapArray( npartTotal, 6 );
		BYTESWAP( flag_cooling );
		BYTESWAP( num_files );
		BYTESWAP( BoxSize );
		BYTESWAP( Omega0 );
		BYTESWAP( OmegaLambda );
		BYTESWAP( HubbleParam );
	}
};

#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3) 			{ if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }



/* This function reads the Dark matter particle information from a general Gadget binary snapshot of type 1 that can be saved in multiple files.

NOTE: This function is a not the best choice to understand how to read input data for users unfamiliar with the Gadget snapshot files.
*/
void readGadgetFile(std::string filename,
                    Read_data<float> *readData,
                    User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string filenameRoot = filename;
    std::string fileName = filename;
    // check to see if there is only one input file or several
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName += "0";
        if ( not bfs::exists(fileName) )
            throwError( "The program could not open the input GADGET snapshot file/files: '" + filenameRoot + "' or '" + fileName + "'. It cannot find the file/files." );
    }
    
    // open the binary file for reading
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );
    
    
    // read the header of the file
    Gadget_header gadgetHeader, tempHeader;
    int gadgetFileType, buffer1, buffer2;   // variable to store gadget file type; variables to read the buffer before and after each gadget data block
	bool swapEndian = false;
    // detect the Gadget file type -> gadget file type 1 or 2
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    if ( buffer1 == 8 )			// gadget file format 2
    {
        gadgetFileType = 2;
        inputFile.seekg( 16, std::ios::beg );
    }
    else if ( buffer1 == 256 )	// gadget file format 1
    {
        gadgetFileType = 1;
        inputFile.seekg( 0, std::ios::beg );
    }
	else						// check for swapped endianness
	{
		BYTESWAP( buffer1 );
		swapEndian = true;
		if ( buffer1 == 8 )			// gadget file format 2
		{
			gadgetFileType = 2;
			inputFile.seekg( 16, std::ios::beg );
		}
		else if ( buffer1 == 256 )	// gadget file format 1
		{
			gadgetFileType = 1;
			inputFile.seekg( 0, std::ios::beg );
		}
		else
			throwError( "Unknown file type for the input Gadget snapshot. Tried Gadget snapshots type 1 and 2 as well as changing endianness, but none worked. Check that you inserted the correct input file." );
		message << "Detected that the input data file has a different endianness than the current system. The program will automatically change endianness for the data!" << MESSAGE::Flush;
	}
    
    
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    inputFile.close();
	SWAP_HEADER_ENDIANNESS( swapEndian, buffer1, buffer2, gadgetHeader ); //swap endianness if that is the case
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET snapshot file. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
    
    // set the box coordinates
    for (size_t i=0; i<NO_DIM; ++i)
    {
        userOptions->boxCoordinates[2*i] = 0.;                       // this is the left extension of the full box
        userOptions->boxCoordinates[2*i+1] = gadgetHeader.BoxSize;   // right extension of the full box
    }
#ifdef WOJTEK
	if ( userOptions->additionalOptions.size()!=0 )	//if inserted an option
		gadgetHeader.num_files = atoi( userOptions->additionalOptions[0].c_str() );
	gadgetHeader.print();
#endif
    
    
    // allocate memory for the particle data and read the data from the binary file
    size_t noParticles = gadgetHeader.npartTotal[1];   //total number of particles - only Dark Matter particles
    Real *positions = readData->position(noParticles); //particle positions
#ifdef VELOCITY
    Real *velocities = readData->velocity(noParticles);//particle velocities
#endif
    Real *weights = readData->weight(noParticles);     //particle weights (weights = particle masses from the snapshot file)
    size_t noReadHalo = 0;                      // the total number of halo particles read after each file
    
    
    // iterate over all the files and read the Dark Matter data in each of them
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        // get the file name
        if ( gadgetHeader.num_files!=1 )
        {
            std::ostringstream buffer;
            buffer << filenameRoot << i;
            fileName = buffer.str();
        }
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;
        
        
        // open the file and read the header
        openInputBinaryFile( inputFile, fileName );
        if ( gadgetFileType==2 ) inputFile.seekg( 16, std::ios::cur );      // jump the block id data if type 2 Gadget file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
		SWAP_HEADER_ENDIANNESS( swapEndian, buffer1, buffer2, tempHeader ); //swap endianness if that is the case
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET file '" + filename + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );
        
        
        // read the position block
        message << "\t reading positions of the particles... " << MESSAGE::Flush;
        if ( gadgetFileType==2 ) inputFile.seekg( 16, std::ios::cur );      // jump the block id data if type 2 Gadget file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        // skip the gas particles
        size_t skipBytes = tempHeader.npart[0]*sizeof(float);
        inputFile.seekg( skipBytes, std::ios::cur );
        // read the halo particles
        size_t readBytes = tempHeader.npart[1]*sizeof(float) * NO_DIM;
        inputFile.read( reinterpret_cast<char *>(&(positions[NO_DIM*noReadHalo])), readBytes );
        // skip the rest of the particles
        skipBytes = (tempHeader.npart[2]+tempHeader.npart[3]+tempHeader.npart[4]+tempHeader.npart[5])*sizeof(float);
        inputFile.seekg( skipBytes, std::ios::cur );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
		SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle position data block in the GADGET file '" + filename + "' did not match. The GADGET snapshot file is corrupt." );
        message << "Done\n";
        
        
#ifdef VELOCITY
        // read the velocity block
        message << "\t reading velocities of the particles... " << MESSAGE::Flush;
        if ( gadgetFileType==2 ) inputFile.seekg( 16, std::ios::cur );      // jump the block id data if type 2 Gadget file
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
         //skip the gas particles
        skipBytes = tempHeader.npart[0]*sizeof(float);
        inputFile.seekg( skipBytes, std::ios::cur );
        // read the halo particles
        readBytes = tempHeader.npart[1]*sizeof(float) * NO_DIM;
        inputFile.read( reinterpret_cast<char *>(&(velocities[3*noReadHalo])), readBytes ); // read the velocities
        // skip the rest of the particles
        skipBytes = (tempHeader.npart[2]+tempHeader.npart[3]+tempHeader.npart[4]+tempHeader.npart[5])*sizeof(float);
        inputFile.seekg( skipBytes, std::ios::cur );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
		SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        if ( buffer1!=buffer2 )
            throwError( "The integers before and after the particle velocity data block in the GADGET file '" + filename + "' did not match. The GADGET snapshot file is corrupt." );
        message << "Done\n";
#else
        // skip the velocity data block if not interested in velocities
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
		SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
#endif
        
        // skiping the identity data block
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
		SWAP_ENDIANNESS( swapEndian, buffer1, buffer2 );
        inputFile.seekg( buffer1+sizeof(int), std::ios::cur );
        
        // read masses if different
        message << "\t reading masses of the particles... " << MESSAGE::Flush;
        if ( tempHeader.mass[1]==0. and tempHeader.npart[1]!=0 )
        {
            if ( gadgetFileType==2 ) inputFile.seekg( 16, std::ios::cur );      // jump the block id data if type 2 Gadget file
            inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
            //skip the gas particles
            if ( tempHeader.mass[0]==0. )
            {
                skipBytes = tempHeader.npart[0]*sizeof(float);
                inputFile.seekg( skipBytes, std::ios::cur );
            }
            readBytes = tempHeader.npart[1]*sizeof(float);
            inputFile.read( reinterpret_cast<char *>(&(weights[noReadHalo])), readBytes ); // read the masses
        }
        else
            for (int i1=0; i1<tempHeader.npart[1]; ++i1)
                weights[noReadHalo+i1] = tempHeader.mass[1];
        message << "Done\n";
        
        
        inputFile.close();
        noReadHalo += tempHeader.npart[1];  //update the number of read DM particles
    }
	
	
	if ( swapEndian )
	{
		ByteSwapArray( positions, NO_DIM*noParticles );
#ifdef VELOCITY
		ByteSwapArray( velocities, NO_DIM*noParticles );
#endif
		if ( tempHeader.mass[1]==0. and tempHeader.npart[1]!=0 )	//if read variable masses
			ByteSwapArray( weights, noParticles );
	}
	
    return;
}


