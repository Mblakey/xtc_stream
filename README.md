# Streaming XTC GROMACS File Decompiler

This is a reverse engineer of the GROMACS XTC format, both from the open-source code and analysis of the binary format. The idea being that we can intercept the decompression process and split into states that allow streaming. Coordinate files can get pretty large (>TiB), so a state saving decompression is a good apprach.

This is a proof of concept, a URL (zenodo) can be inputted which will stream the coordinates to the terminal. Real use case would be streaming the frames into another program, but i'm unaware of what these might be, what the format expected is etc. 


## Build 

Only **curl** is a dependancy here. CMake should find this automatically if system wide. For linux, curl should be packaged already, but can be installed with the distros package manager e.g `sudo apt install curl`. For macos, curl should also be preinstalled but can be manually done with brew e.g `brew install curl`. Once installed, standard cmake build applies:

```
mkdir build
cd build
cmake ..
make
```

## Usage 

```
usage:
  xtc-decompile [options] [<path> | -u <url>]
options:
  -u|--url	stream the xtc file from a valid file url
  -v|--verbose	verbose logging to stderr
```

Common example used for unit testing this was `./xtc-decompile -u "https://zenodo.org/records/14512968/files/example_simulation_6PMB_5000Frames.xtc" > curldecompile.txt` and comparing it to the original GROMACS output from the file on disc. 

## Performance 
