#include <iostream>
#include <cstring>
#include "utilities/OBJLoader.hpp"
#include "utilities/lodepng.h"
#include "rasteriser.hpp"
#include <mpi.h>
#include <string>
#include "chrono"

int main(int argc, char **argv) {
    std::string input("../input/sphere.obj");
    std::string output("../output/sphere.png");
    unsigned int width = 1920;
    unsigned int height = 1080;
    unsigned int depth = 3;


    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    for (int i = 1; i < argc; i++) {
        if (i < argc - 1) {
            if (std::strcmp("-i", argv[i]) == 0) {
                input = argv[i + 1];
            } else if (std::strcmp("-o", argv[i]) == 0) {
                output = argv[i + 1];
                output.erase(output.end() - 4, output.end());
                output += std::to_string(rank) + ".png";
            } else if (std::strcmp("-w", argv[i]) == 0) {
                width = (unsigned int) std::stoul(argv[i + 1]);
            } else if (std::strcmp("-h", argv[i]) == 0) {
                height = (unsigned int) std::stoul(argv[i + 1]);
            } else if (std::strcmp("-d", argv[i]) == 0) {
                depth = (int) std::stoul(argv[i + 1]);
            }
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Loading '" << input << "' file... " << std::endl;

    std::vector<Mesh> meshs = loadWavefront(input, false);

    std::vector<unsigned char> frameBuffer = rasterise(meshs, width, height, depth);

    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time spent rendering the image for process " << rank << " took " << time << " ms" << std::endl;

    MPI_Finalize();

    if (rank == 0) {
        std::cout << "Writing image to '" << output << "'..." << std::endl;
        unsigned error = lodepng::encode(output, frameBuffer, width, height);

        if (error) {
            std::cout << "An error occurred while writing the image file: " << error << ": "
                      << lodepng_error_text(error)
                      << std::endl;
        }
    }

	return 0;
}