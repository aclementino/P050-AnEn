# AnEn - Analog Ensemble Algorithm in C

Computational efficiency improvement of the Analog Ensemble (AnEn) algorithm through C implementation with parallel processing capabilities for multi-core and distributed systems.

## 📝 Project Context

**Institution:** IPB - Polytechnic Institute of Bragança, School of Technology and Management

**Course:** Bachelor's Degree in Computer Engineering (2023/2024)

**Project Theme:** P050 - Improving computational efficiency of the Analog Ensemble algorithm

## 🎯 Project Objectives

The Analog Ensemble (AnEn) algorithm processes large-scale time series, used for weather forecasting or reconstruction of incomplete meteorological data series based on correlated series. This project aims to:

1. Develop a sequential C implementation based on existing R and MATLAB versions
2. Create parallel versions for multi-core systems (using multithreading)
3. Create parallel versions for distributed clusters (using message passing)
4. Ensure scientific validation and publication-ready results

## 📋 Requirements

- GCC 9.0 or higher
- CMake 3.0 or higher
- Make 4.0 or higher
- Build essentials and development libraries

## Dependencies

This project requires the NetCDF-C and GSL (GNU Scientific Library) libraries to be installed before compilation.

## Installation

### Linux (Ubuntu/Debian)

```bash
# Update package repositories
apt-get update

# Install build tools and dependencies
apt-get install -y build-essential git gcc cmake curl m4 libhdf5-dev libxml2 make zlib1g-dev valgrind nano libgsl-dev
apt-get upgrade -y

# Create directory structure for NetCDF-C installation
mkdir -p build/unidata
cd build/unidata

# Clone and install NetCDF-C
git clone https://github.com/Unidata/netcdf-c.git --depth 1
cd netcdf-c
cmake .
make install
make test

# Set library path for NetCDF
export LD_LIBRARY_PATH=/usr/local/lib/

# Add library path to shell configuration (optional, for persistence)
echo 'export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

### Compile Your Project

After NetCDF-C is successfully installed:

```bash
# Make sure the library path is set
export LD_LIBRARY_PATH=/usr/local/lib/

# Build the project
make

# Run the executable
./bin/generic_app

# Clean build files (optional)
make clean

# Rebuild from scratch
make rebuild
```

### Manual Compilation (without Makefile)

If you prefer to compile manually without using the Makefile:

```bash
gcc -O2 -Iinclude -o generic_app src/*.c -lnetcdf -lgsl -lm
```

### Running with Parameters

The application accepts the following parameters:

```bash
./bin/generic_app <threads> <training_period> <netcdf_files>
```

**Parameters:**
- `threads` - Number of threads to use (e.g., 1, 2, 4, 8)
- `training_period` - Training period in years (e.g., 1, 2, 4, 8)
- `netcdf_files` - Path(s) to NetCDF (.nc) files

**Example:**

```bash
./bin/generic_app 2 4 support/nc_data/file1.nc support/nc_data/file2.nc
```

To run multiple tests with different configurations, use the provided test script:

```bash
bash test_script.sh slaptime <max_threads> <num_iterations>
```

This will generate a CSV file with performance metrics in the `test/` directory.

**Note:** The application requires at least 2 NetCDF files to run.

## Troubleshooting

**Error: "command not found: gcc"**
- Install build tools: `sudo apt-get install build-essential`

**Error: "cannot find -lnetcdf"**
- Ensure NetCDF-C installation completed successfully: `make test` should pass
- Verify library path is set: `echo $LD_LIBRARY_PATH`
- Check if libraries exist: `ls -la /usr/local/lib/ | grep netcdf`

**Error: "libnetcdf.so: cannot open shared object file"**
- Set the library path: `export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH`
- For persistent setup, add to `~/.bashrc`

**CMake configuration fails**
- Ensure HDF5 development files are installed: `sudo apt-get install libhdf5-dev`

## Project Structure

```
your-project/
├── src/              # Source files (.c)
├── include/          # Header files (.h)
├── build/            # Build directory (generated)
├── bin/              # Compiled executables (generated)
├── Makefile          # Build configuration
├── README.md         # This file
└── LICENSE           # Project license
```
