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

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen

## at ALBACORE/alba_mod.txt
model = "ALBACORE"

dimension_x = 0 
dimension_y = 0
dimension_z = 0

lon_origin = 0
lat_origin = 0
lon_upper = 0
lat_upper = 0

delta_lon = 0
delta_lat = 0

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url,fname):
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True


def main():

    count =0

    # Set our variable defaults.
    username = ""
    path = ""
    mdir = ""

    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)

    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()

        if (variable == 'model_data_path') :
            path = val + '/' + model
            continue
        if (variable == 'model_dir') :
            mdir = "./"+val
            continue
        if (variable == 'nx') :
            dimension_x = int(val)
            continue
        if (variable == 'ny') :
            dimension_y = int(val)
            continue
        if (variable == 'nz') :
            dimension_z = int(val)
            continue
        if (variable == 'bottom_left_corner_e') :
            lon_origin = float(val)
            continue
        if (variable == 'bottom_left_corner_n') :
            lat_origin = float(val)
            continue
        if (variable == 'top_right_corner_e') :
            lon_upper = float(val)
            continue
        if (variable == 'top_right_corner_n') :
            lat_upper = float(val)
            continue

        continue
    if path == "" :
        print("ERROR: failed to find variables from config file")
        sys.exit(1)

    fp.close()

    delta_lon = (lon_upper - lon_origin )/(dimension_x-1)
    delta_lat = (lat_upper - lat_origin)/(dimension_y-1)
#for ycoord in xrange(dimension_y):
#   yloc = lat_origin + (ycoord * delta_lat)
#  print "yloc - ",yloc
#for xcoord in xrange(dimension_x):
#   xloc = lon_origin + (xcoord * delta_lon)
#   print  "xloc - ",xloc


    try:
        opts, args = getopt.getopt(sys.argv[1:], "", [""])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)

    print("\nDownloading model file\n")

    fname="alba_mod.txt"
    url = path + "/" + fname
    download_urlfile(url,fname)

    # Now we need to go through the data files and put them in the correct
    # format for Albacore. More specifically, we need a Vp.dat, Vs.dat, and Density.dat

    print("\nWriting out ALBACORE data files\n")

    f = open("./alba_mod.txt")

    subprocess.check_call(["mkdir", "-p", mdir])


    f_vp = open(mdir+"/vp.dat", "wb")
    f_vs = open(mdir+"/vs.dat", "wb")
    f_rho = open(mdir+"/rho.dat", "wb")

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


