#!/usr/bin/env python


##
#  Builds the data files in the expected format from alba_mod.txt
#
# from >> lat lon depth(km) vs(km/s) vp(km/s) den(g/cm^2)

# to >>  /** P-wave velocity in meters per second */
#        double vp;
#        /** S-wave velocity in meters per second */
#        double vs;
#        /** Density in g/m^3 */
#        double rho;
#
#   Albacore's rho/density is g/cm^2 ???
 

import getopt
import sys
import subprocess
import struct
import array
import pdb

## at hypocenter  ALBACORE/alba_mod.txt

model = "ALBACORE"
dimension_x = 27
dimension_y = 8
dimension_z = 101

lon_origin = -124.6472
lat_origin = 32.7

lon_upper = -116.84720
lat_upper = 34.8

delta_lon = (lon_upper - lon_origin )/(dimension_x-1)
delta_lat = (lat_upper - lat_origin)/(dimension_y-1)

for ycoord in xrange(dimension_y):
   yloc = lat_origin + (ycoord * delta_lat)
#  print "yloc - ",yloc
for xcoord in xrange(dimension_x):
   xloc = lon_origin + (xcoord * delta_lon)
#   print  "xloc - ",xloc

def usage():
    print("\n./make_data_files.py -u [uid]\n\n")
    print("-u - username to use to do the dataset retrieval.\n")
    sys.exit(0)

def main():

    count =0

    # Set our variable defaults.
    username = ""
    path = "/var/www/html/research/ucvmc/" + model 

    try:
        opts, args = getopt.getopt(sys.argv[1:], "u:", ["user="])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ("-u", "--user"):
            username = str(a) + "@"

    print("\nDownloading model file\n")

##
##    subprocess.check_call(["scp", username +
##                           "hypocenter.usc.edu:" + path + "/alba_mod.txt",
##                           "."])
##

    # Now we need to go through the data files and put them in the correct
    # format for Albacore. More specifically, we need a Vp.dat, Vs.dat, and Density.dat

    print("\nWriting out ALBACORE data files\n")

    f = open("./alba_mod.txt")

    f_vp = open("./alba/vp.dat", "wb")
    f_vs = open("./alba/vs.dat", "wb")
    f_rho = open("./alba/rho.dat", "wb")

    vp_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))
    vs_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))
    rho_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))

    first = True
    for line in f:
        arr = line.split()
        if first :  ## skip first line
           first = False
           continue
        lat_v = float(arr[0])
        lon_v = float(arr[1])
        depth_v = float(arr[2])
        vs = float(arr[3]) * 1000 
        vp = float(arr[4]) * 1000
        rho = float(arr[5]) * 1000

        y_pos = int(round((lat_v - lat_origin) / delta_lat))
        x_pos = int(round((lon_v - lon_origin) / delta_lon))
        z_pos = int(depth_v)

        count = count + 1
        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos

        vp_arr[loc] = vp
        vs_arr[loc] = vs
        rho_arr[loc] = rho
        if( loc < 10) :
           print line
           print "XXX loc", loc, ":vp", vp," vs",vs," rho",rho
           print " ==> ", x_pos," ",y_pos," ",z_pos," >> ", lon_v, " ",lat_v, " ",depth_v 

    vp_arr.tofile(f_vp)
    vs_arr.tofile(f_vs)
    rho_arr.tofile(f_rho)

    f.close()
    f_vp.close()
    f_vs.close()
    f_rho.close()

    print("\nDone!")

if __name__ == "__main__":
    main()


