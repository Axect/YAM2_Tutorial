from netCDF4 import Dataset
import numpy as np

g = np.array([1,-1,-1,-1])

def np_mass(p4, ax=1):
    return np.sqrt(np.abs(np.sum(g*(p4**2),axis=ax)))

toy_data = np.load("data/toy_array_new.npz")

pa1 = toy_data['b1']
pa2 = toy_data['b2']

pb1 = toy_data['l2']
pb2 = toy_data['l1']

nu_1 = toy_data['nu1']
nu_2 = toy_data['nu2']

print(pa1[0:5,0])

print(pa1.shape)

ncfile = Dataset("data/toy_array_new.nc", 'w', format='NETCDF4')
ncfile.createDimension('np', pa1.shape[0])

pa1_var_E = ncfile.createVariable('pa1_E', 'float64', ('np')); pa1_var_E[:] = pa1[:,0]
pa1_var_x = ncfile.createVariable('pa1_x', 'float64', ('np')); pa1_var_x[:] = pa1[:,1]
pa1_var_y = ncfile.createVariable('pa1_y', 'float64', ('np')); pa1_var_y[:] = pa1[:,2]
pa1_var_z = ncfile.createVariable('pa1_z', 'float64', ('np')); pa1_var_z[:] = pa1[:,3]

pa2_var_E = ncfile.createVariable('pa2_E', 'float64', ('np')); pa2_var_E[:] = pa2[:,0]
pa2_var_x = ncfile.createVariable('pa2_x', 'float64', ('np')); pa2_var_x[:] = pa2[:,1]
pa2_var_y = ncfile.createVariable('pa2_y', 'float64', ('np')); pa2_var_y[:] = pa2[:,2]
pa2_var_z = ncfile.createVariable('pa2_z', 'float64', ('np')); pa2_var_z[:] = pa2[:,3]

pb1_var_E = ncfile.createVariable('pb1_E', 'float64', ('np')); pb1_var_E[:] = pb1[:,0]
pb1_var_x = ncfile.createVariable('pb1_x', 'float64', ('np')); pb1_var_x[:] = pb1[:,1]
pb1_var_y = ncfile.createVariable('pb1_y', 'float64', ('np')); pb1_var_y[:] = pb1[:,2]
pb1_var_z = ncfile.createVariable('pb1_z', 'float64', ('np')); pb1_var_z[:] = pb1[:,3]

pb2_var_E = ncfile.createVariable('pb2_E', 'float64', ('np')); pb2_var_E[:] = pb2[:,0]
pb2_var_x = ncfile.createVariable('pb2_x', 'float64', ('np')); pb2_var_x[:] = pb2[:,1]
pb2_var_y = ncfile.createVariable('pb2_y', 'float64', ('np')); pb2_var_y[:] = pb2[:,2]
pb2_var_z = ncfile.createVariable('pb2_z', 'float64', ('np')); pb2_var_z[:] = pb2[:,3]

nu_1_var_E = ncfile.createVariable('nu1_E', 'float64', ('np')); nu_1_var_E[:] = nu_1[:,0]
nu_1_var_x = ncfile.createVariable('nu1_x', 'float64', ('np')); nu_1_var_x[:] = nu_1[:,1]
nu_1_var_y = ncfile.createVariable('nu1_y', 'float64', ('np')); nu_1_var_y[:] = nu_1[:,2]
nu_1_var_z = ncfile.createVariable('nu1_z', 'float64', ('np')); nu_1_var_z[:] = nu_1[:,3]

nu_2_var_E = ncfile.createVariable('nu2_E', 'float64', ('np')); nu_2_var_E[:] = nu_2[:,0]
nu_2_var_x = ncfile.createVariable('nu2_x', 'float64', ('np')); nu_2_var_x[:] = nu_2[:,1]
nu_2_var_y = ncfile.createVariable('nu2_y', 'float64', ('np')); nu_2_var_y[:] = nu_2[:,2]
nu_2_var_z = ncfile.createVariable('nu2_z', 'float64', ('np')); nu_2_var_z[:] = nu_2[:,3]

ncfile.close();