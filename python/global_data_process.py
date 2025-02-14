"""
Created on Sun Jun 4 15:45:29 2023

Author: Wei-Ting Hung
"""

import os
import subprocess
from datetime import datetime, timezone

import numpy as np
from netCDF4 import Dataset
from pysolar.solar import get_altitude
from scipy.interpolate import griddata

"""User Options"""
path = "/scratch/pcampbe8/canopy-app/input"  # work directory
year = 2022  # year
month = 7  # month
day = 1  # day
houri = 12  # gfs initialization hour in UTC (caution: currently GFS input files are initialized at 12 UTC only -- do not change)
hour = 0  # gfs forecast hour (0-24) in UTC
ref_lev = 10  # reference height above the canopy (m)
frp_src = 1  # frp data source for WAF (0: local source; 1: check local source first, switch to climatological file if no available data; 2: 12 month climatology; 3: all ones when ifcanwaf=.FALSE.)


# ---------------------------------ATTENTION---------------------------------#
# If local FRP is used (frp_src=0,1), archived GBBEPx files since 2020 are   #
# available for GMU HOPPER users. For users outside GMU, parameter "f_frp"   #
# and function "read_frp_local" need to be modified accordingly.             #
# Function read_frp_local is designed for reading GBBEPx by default.         #
#                                                                            #
# Recent GBBEPx v3 files (~ < 3 months) are available for download at:       #
# https://www.ospo.noaa.gov/Products/land/gbbepx/                            #
# GBBEPx_all01GRID.emissions_v003_'+YY+MM+DD+'.nc'                           #
#                                                                            #
# Archived GBBEPx files since 2020 are available for GMU HOPPER users.       #
#                                                                            #
# 12 month climatological FRP will not be respective of actual conditions.   #
# Only use it when users do not need the limited WAF application for fires.  #
# ---------------------------------------------------------------------------#


starttime = datetime.now()
print("------------------------------------")
print("---- Global input pre-process start!", starttime.strftime("%Y/%m/%d %H:%M:%S"))
print("------------------------------------")


"""Settings"""
# date/time
YY = str(year)
MM = "%02d" % month
DD = "%02d" % day
HHI = "%02d" % houri
HH = "%02d" % hour


# domain
lat_lim = [-90, 90]
lon_lim = [0, 360]


# input/output files
f_met = (
    path + "/gfs.t" + HHI + "z." + YY + MM + DD + ".sfcf0" + HH + ".nc"
)  # gfs met file
f_can = path + "/gfs.canopy.t" + HHI + "z." + YY + MM + "01.sfcf000.nc"  # canopy file
f_output = (
    path + "/gfs.canopy.t" + HHI + "z." + YY + MM + DD + ".global.f0" + HH + ".nc"
)  # output file

if (frp_src == 0) or (frp_src == 1):  # local frp file
    if int(YY + MM + DD) <= 20230510:  # version 3
        f_frp = (
            "/groups/ESS/yli74/data/GBBEPx/ORI/GBBEPx_all01GRID.emissions_v003_"
            + YY
            + MM
            + DD
            + ".nc"
        )
    else:  # version 4
        f_frp = (
            "/groups/ESS/yli74/data/GBBEPx/ORI/GBBEPx-all01GRID_v4r0_blend_"
            + YY
            + MM
            + DD
            + ".nc "
        )
elif frp_src == 2:  # climatological frp
    f_frp = path + "/gfs.canopy.t" + HHI + "z." + YY + MM + "01.sfcf000.nc"


# required variables
loclist = ["grid_xt", "lon", "grid_yt", "lat", "time"]
metlist = [
    "ugrd10m",
    "vgrd10m",
    "fricv",
    "sfcr",
    "vtype",
    "sotyp",
    "pressfc",
    "dswrf",
    "shtfl",
    "tmpsfc",
    "tmp2m",
    "spfh2m",
    "hpbl",
    "prate_ave",
]
canlist = ["lai", "clu", "ffrac", "fh", "mol", "csz", "frp", "href"]


# constants
fill_value = 9.99e20  # fill value
den = 1.18  # air density (kg/m3)
Cp = 1004  # specific heat capacity of Air at 25C (J/kg/K)
K = 0.4  # Von Karman constant
g = 9.8  # gravitational acceleration (m/s2)


# functions
def read_varatt(var):
    attname = var.ncattrs()
    att = [var.getncattr(X) for X in attname]
    return attname, att


def write_varatt(var, attname, att):
    for X in np.arange(len(attname)):
        if attname[X] == "_FillValue":
            continue
        elif attname[X] == "missing_value":
            value = np.float32(att[X])
            value = np.round(value / 1e15) * 1e15
            var.setncattr(attname[X], value)
        else:
            var.setncattr(attname[X], att[X])


def mapping(xgrid, ygrid, data, xdata, ydata, map_method, fvalue):
    output = griddata(
        (xdata, ydata), data, (xgrid, ygrid), method=map_method, fill_value=fvalue
    )
    return output


def read_gfs_climatology(filename, lat, lon, varname):
    readin = Dataset(filename)

    # map to met grids
    yt = readin["lat"][:]
    xt = readin["lon"][:]
    data = np.squeeze(readin[varname][0, :, :])

    DATA = mapping(lat, lon, data.flatten(), yt.flatten(), xt.flatten(), "linear", np.nan)
    DATA[np.isnan(DATA)] = 0
    DATA[DATA < 0] = 0
    return DATA


def read_frp_local(filename, lat, lon, fill_value):
    readin = Dataset(filename)

    # map to met grids
    xt, yt = np.meshgrid(readin["Longitude"][:], readin["Latitude"][:])
    xt[xt < 0] = xt[xt < 0] + 360
    data = np.squeeze(readin["MeanFRP"][:])

    DATA = mapping(
        lat, lon, data.flatten(), yt.flatten(), xt.flatten(), "nearest", fill_value
    )
    return DATA


"""Data Download"""
"""Download from servers if required files do not exist."""
print("---- Checking required files...")
print("------------------------------------")

# met file
if os.path.isfile(f_met) is True:
    print("---- Met file found!")
else:
    print("---- Cannot find met file. Downloading from AWS...")
    subprocess.run(
        [
            "wget",
            "--no-check-certificate",
            "--no-proxy",
            "-O",
            path + "/gfs.t" + HHI + "z." + YY + MM + DD + ".sfcf0" + HH + ".nc",
            "https://nacc-in-the-cloud.s3.amazonaws.com/inputs/"
            + YY
            + MM
            + DD
            + "/gfs.t"
            + HHI
            + "z.sfcf0"
            + HH
            + ".nc",
        ]
    )
    if os.path.isfile(f_met) is True:
        os.chmod(f_met, 0o0755)
        print("---- Download complete!")
    else:
        print("---- No available met data. Terminated!")
        exit()

# can file
if os.path.isfile(f_can) is True:
    print("---- Canopy file found!")
else:
    print("---- Cannot find canopy file. Downloading from AWS...")
    subprocess.run(
        [
            "wget",
            "--no-check-certificate",
            "--no-proxy",
            "-P",
            path,
            "https://nacc-in-the-cloud.s3.amazonaws.com/inputs/geo-files/gfs.canopy.t"
            + HHI
            + "z."
            + YY
            + MM
            + "01.sfcf000.nc",
        ]
    )
    if os.path.isfile(f_can) is True:
        os.chmod(f_can, 0o0755)
        print("---- Download complete!")
    else:
        print("---- No available canopy data. Terminated!")
        exit()

# frp file
if frp_src == 0:  # local source
    if os.path.isfile(f_frp) is True:
        os.system("cp " + f_frp + " " + path)
        if int(YY + MM + DD) <= 20230510:
            f_frp = path + "/GBBEPx_all01GRID.emissions_v003_" + YY + MM + DD + ".nc"
        else:
            f_frp = path + "/GBBEPx-all01GRID_v4r0_blend_" + YY + MM + DD + ".nc"
        os.chmod(f_frp, 0o0755)
        print("---- FRP file found!")
    else:
        print("---- No available FRP file. Terminated!")
        exit()

if frp_src == 1:  # local source
    if os.path.isfile(f_frp) is True:
        os.system("cp " + f_frp + " " + path)
        if int(YY + MM + DD) <= 20230510:
            f_frp = path + "/GBBEPx_all01GRID.emissions_v003_" + YY + MM + DD + ".nc"
        else:
            f_frp = path + "/GBBEPx-all01GRID_v4r0_blend_" + YY + MM + DD + ".nc"
        os.chmod(f_frp, 0o0755)
        print("---- FRP file found!")
    else:
        print("---- No available FRP file. Switch to Climatology FRP...")
        frp_src = 2
        f_frp = path + "/gfs.canopy.t" + HHI + "z." + YY + MM + "01.sfcf000.nc"

if frp_src == 2:  # 12 month climatology frp
    if os.path.isfile(f_frp) is True:
        print("---- FRP file found!")
    else:
        print("---- Canot find FRP file. Downloading from AWS...")
        subprocess.run(
            [
                "wget",
                "--no-check-certificate",
                "--no-proxy",
                "-P",
                path,
                "https://nacc-in-the-cloud.s3.amazonaws.com/inputs/geo-files/gfs.canopy.t"
                + HHI
                + "z."
                + YY
                + MM
                + "01.sfcf000.nc",
            ]
        )
        if os.path.isfile(f_met) is True:
            os.chmod(f_frp, 0o0755)
            print("---- Download complete!")
        else:
            print("---- No available FRP data. Terminated!")
            exit()

    print("-----------!!!WARNING!!!-------------")
    print("---!!!Climatological FRP is used!!!--")


os.system("cp " + f_met + " " + f_output)  # copy gfs met file for appending


"""Reading dimensions"""
print("------------------------------------")
print("---- Checking variable dimensions...")
print("------------------------------------")
readin = Dataset(f_met)
grid_yt = readin["grid_yt"][:]
grid_xt = readin["grid_xt"][:]
lat = readin["lat"][:]
lon = readin["lon"][:]
time = readin["time"][:]


# dimension sizes
ntime = len(time)
nlat = len(grid_yt)
nlon = len(grid_xt)


# var check
print("time", time.shape)
print("grid_yt", grid_yt.shape)
print("grid_xt", grid_xt.shape)
print("lat", lat.shape)
print("lon", lon.shape)


"""Adding canvar"""
print("------------------------------------")
print("---- Generating canopy variables...")
print("------------------------------------")

for i in np.arange(len(canlist)):
    varname = canlist[i]

    print("---- " + varname + " processing...")

    if varname == "lai":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Leaf area index", "m^2/m^2", fill_value]
        DATA = read_gfs_climatology(f_can, lat, lon, "lai")

    elif varname == "clu":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Canopy clumping index", "none", fill_value]
        DATA = read_gfs_climatology(f_can, lat, lon, "clu")

    elif varname == "ffrac":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Forest fraction of grid cell", "none", fill_value]
        DATA = read_gfs_climatology(f_can, lat, lon, "ffrac")

    elif varname == "fh":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Canopy height above the surface", "m", fill_value]
        DATA = read_gfs_climatology(f_can, lat, lon, "fh")

    elif varname == "mol":
        # Reference:
        # Essa 1999, ESTIMATION OF MONIN-OBUKHOV LENGTH USING RICHARDSON AND BULK RICHARDSON NUMBER
        # https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/118/37118528.pdf
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Monin-Obukhov length", "m", fill_value]

        readin = Dataset(f_met)
        t2m = np.squeeze(readin["tmp2m"][:])
        fricv = np.squeeze(readin["fricv"][:])
        shtfl = np.squeeze(readin["shtfl"][:])

        DATA = (-1 * den * Cp * t2m * (fricv**3)) / (K * g * shtfl)
        DATA[DATA > 500] = 500
        DATA[DATA < -500] = -500

        del [readin, t2m, fricv, shtfl]

    elif varname == "csz":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Cosine of solar zenith angle", "none", fill_value]

        time_conv = datetime(
            int(YY), int(MM), int(DD), int(HH), 0, 0, 0, tzinfo=timezone.utc
        )
        sza = 90 - get_altitude(lat, lon, time_conv)
        DATA = np.cos(sza * 0.0174532925)  # degree to radian

        del [time_conv, sza]

    elif varname == "frp":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Mean fire radiative power", "MW", fill_value]

        if frp_src == 2:
            DATA = read_gfs_climatology(f_can, "frp")
        elif frp_src == 3:
            DATA = np.empty(lat.shape)
            DATA[:] = 1
        else:
            DATA = read_frp_local(f_frp, lat, lon, fill_value)

    elif varname == "href":
        ATTNAME = ["long_name", "units", "missing_value"]
        ATT = ["Reference height above the canopy", "m", fill_value]
        DATA = np.empty([nlat, nlon])
        DATA[:] = ref_lev

    # var check
    print("Dimension/Attributes:")
    print(DATA.shape)
    print(ATTNAME)
    print(ATT)

    # adding to output file
    output = Dataset(f_output, "a")

    var = output.createVariable(
        varname, "float", ("time", "grid_yt", "grid_xt"), fill_value=fill_value
    )
    write_varatt(var, ATTNAME, ATT)
    var[:] = DATA

    output.close()

    print("---- " + varname + " complete!")

    del [output, var]
    del [varname, DATA, ATTNAME, ATT]


endtime = datetime.now()


print("------------------------------------")
print("---- Global input pre-process complete!", endtime.strftime("%Y/%m/%d %H:%M:%S"))
print("---- Process time:", str(endtime - starttime))
print("------------------------------------")
