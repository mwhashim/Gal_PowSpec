#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

void readBinaryFile(std::string filename,
                    Read_data<float> *readData)
{
    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, filename );
    
    
    // read the number of particles and the box coordinates from the binary file
    int noParticles;
    float boxCoordinates[2*NO_DIM];
    inputFile.read( reinterpret_cast<char *>(&noParticles), sizeof(noParticles) );
    inputFile.read( reinterpret_cast<char *>(boxCoordinates), sizeof(float) );
    for (size_t i=0; i<2*NO_DIM; ++i)
        userOptions->boxCoordinates[i] = boxCoordinates[i];
    
    
    // reserve memory for the input data
    Real *positions = readData->position(noParticles);  //particle positions
    Real *weights = readData->weight(noParticles);      //particle weights (e.g. weights = particle masses)
    Real *velocities = readData->velocity(noParticles); //particle velocities
    
    
    // read the rest of the input data: positions, weights and velocities
    size_t dataSize = noParticles * sizeof(float) * NO_DIM;    // number of data bytes that store the particle positions (3*4 bytes per particle)
    inputFile.write( reinterpret_cast<char const *>(positions), dataSize );
    
    dataSize = noParticles * sizeof(float);    // number of data bytes that store the particle weights (1*4 bytes per particle)
    inputFile.write( reinterpret_cast<char const *>(positions), dataSize );
    
    dataSize = noParticles * sizeof(float) * NO_DIM;    // number of data bytes that store the particle velocities (3*4 bytes per particle)

		inputFile.write( reinterpret_cast<char const *>(positions), dataSize );
    
    inputFile.close();
    
}


int main(){

	return 1;
}


