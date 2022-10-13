import numpy as np
import matplotlib.pyplot as plt
import h5py

g = np.array([1,-1,-1,-1])

def np_mass(p4, ax=1):
    return np.sqrt(np.abs(np.sum(g*(p4**2),axis=ax)))

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

with h5py.File('data/m2cc.h5', 'r') as df:
    print("Keys: %s" % list(df.keys()))
    m2cc = np.array(list(df['m2cc']))
    k1 = np.array(list(df['k1']))
    k2 = np.array(list(df['k2']))
    pa1 = np.array(list(df['pa1']))
    pa2 = np.array(list(df['pa2']))
    pb1 = np.array(list(df['pb1']))
    pb2 = np.array(list(df['pb2']))

k1 = k1.reshape(4, m2cc.shape[0]).T
k2 = k2.reshape(4, m2cc.shape[0]).T
pa1 = pa1.reshape(4, m2cc.shape[0]).T
pa2 = pa2.reshape(4, m2cc.shape[0]).T
pb1 = pb1.reshape(4, m2cc.shape[0]).T
pb2 = pb2.reshape(4, m2cc.shape[0]).T

mc1 = k1[:,0]**2 - k1[:,1]**2 - k1[:,2]**2 - k1[:,3]**2
mc2 = k2[:,0]**2 - k2[:,1]**2 - k2[:,2]**2 - k2[:,3]**2

# Prepare Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$M_2$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(m2cc, bins=2000, label='YAM2')

# Other options
plt.xlim(left=600, right=1200)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("m2cc_hist.png", dpi=300)

# ==============================================================================
# m_B
# ==============================================================================
pB1 = pb1 + k1
pB2 = pb2 + k2

mB1 = np_mass(pB1)
mB2 = np_mass(pB2)

# Prepare Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$M_{B_1}$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(mB1, bins=2000, label='YAM2')

# Other options
plt.xlim(left=600, right=1000)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("mB1_hist.png", dpi=300)

# Prepare Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$M_{B_2}$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(mB2, bins=2000, label='YAM2')

# Other options
plt.xlim(left=600, right=1000)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("mB2_hist.png", dpi=300)

# ==============================================================================
# m_A
# ==============================================================================
pA1 = pa1 + pB1
pA2 = pa2 + pB2

mA1 = np_mass(pA1)
mA2 = np_mass(pA2)

# Prepare Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$M_{A_1}$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(mA1, bins=2000, label='YAM2')

# Other options
plt.xlim(left=600, right=1200)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("mA1_hist.png", dpi=300)

# Prepare Plot
plt.figure(figsize=(10,6), dpi=300)
plt.title(r"$M_{A_2}$", fontsize=16)
plt.xlabel(r'Mass [GeV]', fontsize=14)
plt.ylabel(r'Number', fontsize=14)

# Plot with Legends
plt.hist(mA2, bins=2000, label='YAM2')

# Other options
plt.xlim(left=600, right=1200)
plt.legend(fontsize=12)
plt.grid()
plt.savefig("mA2_hist.png", dpi=300)
