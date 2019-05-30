from netCDF4 import Dataset
import sys
import os

if len(sys.argv) <= 1 :
    print("Error: Missing filename argument")
    print("Usage: python convert_to_netcdf4.py filename")

fname_in  = sys.argv[1]
fname_out = sys.argv[1].replace(".nc","_nc4.nc")

src = Dataset(fname_in ,"r")
dst = Dataset(fname_out,"w",format="NETCDF4")

# copy dimensions
for name, dimension in src.dimensions.items():
    dst.createDimension(
        name, (len(dimension) if not dimension.isunlimited() else None))
# copy all file data except for the excluded
for name, variable in src.variables.items():
    x = dst.createVariable(name, variable.datatype, variable.dimensions)
    dst[name][:] = src[name][:]

src.close()
dst.close()
