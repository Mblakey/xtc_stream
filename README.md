# Streaming XTC GROMACS File Decompiler

This is a reverse engineer of the GROMACS XTC format, both from the open-source code and analysis of the binary format. The idea being that we can intercept the decompression process and split into states that allow streaming. 

Josh mentioned that coordinate files can get pretty chunky (>TiB) and that it would be a nice idea to stream these into other programs so that the memory footprint doesn't exceed a certain size. 

This is a proof of concept, a URL (zenodo) can be inputted which will stream the coordinates to the terminal. Real use case would be streaming the frames into another program, but i'm unaware of what these might be, what the format expected is etc. 

For any other C developers looking at my code, the switch state and fallthrough label model was used as it plays really nicely with gdb when reverse engineering the format on the fly. It's not built for performance. I will touch base with what we need before optimising this (even the original C code wasn't great.. tut tut GROMACS). 

Low hanging optimisations would be spooling the curl buffer, removing the switch labels for functions or goto jmp statements, and then using multithreading curl functions to increase the download speed (this will likely bottleneck the process). 

If needed, a socket to a IPv4|6 address for a LAN Gigabit connection would give the best performance if trying to stream from a workstation to a laptop. 

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

Common example I used for unit testing this was `./xtc-decompile -u "https://zenodo.org/records/14512968/files/example_simulation_6PMB_5000Frames.xtc" > curldecompile.txt` and comparing it to the original GROMACS output from the file on disc. 


Anyway.. hope this is a useful proof of concept! have fun Josh. 

Michael Blakey 
Senior Software Engineer (Systems & Performance)
NextMove Software

