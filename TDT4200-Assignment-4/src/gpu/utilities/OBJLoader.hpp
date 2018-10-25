#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include "geometry.hpp"
#include "cuda_runtime.h"

std::vector<GPUMesh> loadWavefrontGPU(std::string const srcFile, bool quiet = true);