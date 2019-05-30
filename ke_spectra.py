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

u = ncid.variables['u'][nt-2,nz/2,:,:]
v = ncid.variables['v'][nt-2,nz/2,:,:]
w = ncid.variables['w'][nt-2,nz/2,:,:]

for j in range(ny) :
  u[j,:] = u[j,:] - sum(u[j,:])/nx
  v[j,:] = v[j,:] - sum(v[j,:])/nx
  w[j,:] = w[j,:] - sum(w[j,:])/nx

print( np.amax(u) )
print( np.amax(v) )
print( np.amax(w) )

ke = (u*u+v*v+w*w)/2

#Compute a y-average of the x-direction FFTs 
sp = np.fft.rfft(ke[0,:])[0:int(nx/2)]
sp[:] = 0
for j in range(ke.shape[1]) :
    sp = sp + np.fft.rfft(ke[j,:])[0:int(nx/2)] / ke.shape[1];
spd = sp
for i in range(sp.shape[0]) :
    spd[i] = np.real(sp[i])*np.real(sp[i]) + np.imag(sp[i])*np.imag(sp[i])
spd[1:] = spd[1:]*2

#Compute the wavenumbers associated with the FFT
n_value = np.fft.fftfreq( nx , (1.0 / (nx) ) )
wn = np.array([0. for i in range(spd.shape[0])])
for i in range(int(nx/2)) :
    wn[i] = 2*np.pi*n_value[i]/ncid.variables['x'][nx-1]

spd.dump('spd.npy')
wn.dump('wn.npy')

slope = 1e-1*np.power(wn,-5./3.)

plt.loglog(wn,spd)
plt.loglog(wn,slope)
plt.xlabel('Wavenumber')
plt.ylabel('Kinetic Energy Spectral Power Density')
wn6dx = 2*np.pi/(6*20000./nx)
wn4dx = 2*np.pi/(4*20000./nx)
wn2dx = 2*np.pi/(2*20000./nx)
plt.axvline(x=wn2dx,linestyle="--")
plt.axvline(x=wn4dx,linestyle="--")
# plt.axvline(x=wn6dx,linestyle="--")
plt.show()


# plt.savefig('collision_ke_spectra.eps', bbox_inches='tight')
# plt.close()
