# rad2ply
Glacier radiometry data to Stanford PLY converter

This tool converts radiometry data into simple geometry. Return power is coded into vertex colors.

## Prerequisites

Requires `pickle` (built-in), `NumPy` and `osgeo` modules.

## Usage

Right now, input and output files and various parameters are hard coded in `makeply.py` . Make changes as desired and then simply run

    ./makeply.py
   
If all goes well, you will have a new .ply file that you can load into, say, [Blender](https://www.blender.org).
