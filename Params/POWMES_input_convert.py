from __future__ import division
from numpy import *

x, y, z, lg10M, v_x, v_y, v_z  = loadtxt("/users/mhashim/lustre/mhashim/OxProj/Galaxy_Catalog/tmpcat_761_laigle.dat", unpack=True)
h = 0.704
x_new = (x*h + 50.0)/100.0; y_new = (y*h + 50.0)/100.0; z_new = (z*h + 50.0)/100.0

f = open("/users/mhashim/lustre/mhashim/OxProj/Galaxy_Catalog/tmpcat_761_laigle_inputPOWMES.dat", "w")
f.writelines("%s\n" %(126361))
f.writelines("%f	%f	%f	%f\n" %(x_new[i], y_new[i], z_new[i], 10**(lg10M[i] - 10.0)*h) for i in range(len(x)))
f.close()
