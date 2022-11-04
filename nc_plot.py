from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# Open the netCDF file
nc = Dataset('./data/m2ccb.nc')
var = nc.variables
m2ccb = np.array(var['m2'][:])

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Prepare Plot
plt.figure(figsize=(10,6), dpi=150)
plt.title(r"$M_2$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(m2ccb, bins=2000, label=r'$M_2^{CC}(b)$', histtype='step', linewidth=2)

# Other options
plt.xlim(left=600, right=1000)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("m2ccb_hist.png", dpi=300)

# Open the netCDF file
nc = Dataset('./data/ttbar_m2ccb.nc')
var = nc.variables
m2ccb = np.array(var['m2'][:])

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Prepare Plot
plt.figure(figsize=(10,6), dpi=150)
plt.title(r"$M_2 (t\bar{t})$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(m2ccb, bins=2000, label=r'$M_2^{CC}(b)$', histtype='step', linewidth=2)

# Other options
plt.xlim(left=0, right=150)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("ttbar_m2ccb_hist.png", dpi=300)
