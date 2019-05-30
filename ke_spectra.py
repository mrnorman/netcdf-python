import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

def running_mean(x, N):
    ret = x
    h = int(np.floor(N/2))
    nx = x.shape[0]
    for i in range(nx) :
        if (i-h<0) :
            ret[i] = np.average(x[0:i+1])
        elif (i-h+N>nx-1) :
            ret[i] = np.average(x[i:nx])
        else :
            ret[i] = np.average(x[i-h:i-h+N])
    return ret


###################################################################
## COLLISION
###################################################################
ncid = nc.Dataset("output.nc", "r")
nx = ncid.dimensions['x'].size
ny = ncid.dimensions['y'].size
nz = ncid.dimensions['z'].size
nt = ncid.dimensions['t'].size
x,y = np.meshgrid(ncid.variables['x'][:]/1000.,ncid.variables['y'][:]/1000.)
t = ncid.variables['t'][:]

print(nx)
print(ny)
print(nz)
print(nt)

u = ncid.variables['u'][nt-1,3*nz/4,:,:]
v = ncid.variables['v'][nt-1,3*nz/4,:,:]
w = ncid.variables['w'][nt-1,3*nz/4,:,:]

u = u - sum(u)/nx/ny
v = v - sum(u)/nx/ny
w = w - sum(u)/nx/ny

ke = (u*u+v*v+w*w)/2

#Compute a z-average of the x-direction FFTs 
sp = np.fft.rfft(ke[0,:])[0:int(nx/2)]
sp[:] = 0
for k in range(ke.shape[0]) :
    sp = sp + np.fft.rfft(ke[k,:])[0:int(nx/2)]
spd = sp
for i in range(sp.shape[0]) :
    spd[i] = np.real(sp[i])*np.real(sp[i]) + np.imag(sp[i])*np.imag(sp[i])
spd[1:] = spd[1:]*2

#Compute the wavenumbers associated with the FFT
n_value = np.fft.fftfreq( nx , (1.0 / (nx) ) )
wn = np.array([0. for i in range(spd.shape[0])])
for i in range(int(nx/2)) :
    wn[i] = 2*np.pi*n_value[i]/ncid.variables['x'][nx-1]

plt.loglog(wn,running_mean(spd,3))
plt.xlabel('Wavenumber')
plt.ylabel('Kinetic Energy Spectral Power Density')
plt.savefig('collision_ke_spectra.eps', bbox_inches='tight')
plt.close()
