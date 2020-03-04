#!/usr/bin/env python

##
#  Retrieve the latlons and retrieve the material properties
#

import getopt
import sys
import subprocess
import struct
import numpy as np
import pdb

dimension_x = 27
dimension_y = 8
dimension_z = 101

lon_origin = -124.6472
lat_origin = 32.7

lon_upper = -116.84720
lat_upper = 34.8

delta_lon = (lon_upper - lon_origin )/(dimension_x-1)
delta_lat = (lat_upper - lat_origin)/(dimension_y-1)

def usage():
    print("\n./query_data_files.py\n\n")
    sys.exit(0)

def myrange(start, end, step):
    while start <= end:
        yield start
        start += step

def main():

    count =0

    f_lonlat = open("./lonlat.txt")

    f_vp = open("./alba/vp.dat")
    f_vs = open("./alba/vs.dat")
    f_density = open("./alba/density.dat")

    vp_arr = np.fromfile(f_vp, dtype=np.float32)
    vs_arr = np.fromfile(f_vs, dtype=np.float32)
    density_arr = np.fromfile(f_density, dtype=np.float32)

    f_vp.close()
    f_vs.close()
    f_density.close()

    lon_start = lon_origin
    lat_start = lat_origin
    depth_start = 0.0;

    for line in f_lonlat:
        arr = line.split()
        lon_v = float(arr[0])
        lat_v = float(arr[1])
        depth_v = float(arr[2])

        y_pos = int(round((lat_v - lat_origin) / delta_lat))
        x_pos = int(round((lon_v - lon_origin) / delta_lon))
        z_pos = int(depth_v)

        offset=z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vp=vp_arr[offset];
        vs=vs_arr[offset];
        density=density_arr[offset];

        print x_pos," ",y_pos," ",z_pos," >> ", lon_v, " ",lat_v, " ", float(depth_v) , "-->", vp," ", vs," ", density

    f_lonlat.close()
    print("\nDone!")

if __name__ == "__main__":
    main()


