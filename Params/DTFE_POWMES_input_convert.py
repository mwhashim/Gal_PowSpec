from __future__ import division
from numpy import *

x, y, z, delta_M  = loadtxt("/users/mhashim/lustre/mhashim/OxProj/OutPut/tmpcat_761_laigle_txt_32.a_den", unpack=True)
x_new = (x + 50.0)/100.0; y_new = (y + 50.0)/100.0; z_new = (z + 50.0)/100.0
f = open("/users/mhashim/lustre/mhashim/OxProj/Galaxy_Catalog/tmpcat_761_laigle_inputPOWMES_DTFE.dat", "w")
f.writelines("%s\n" %(32768))
f.writelines("%f	%f	%f	%f\n" %(x_new[i], y_new[i], z_new[i], delta_M[i]) for i in range(len(x)))
f.close()
