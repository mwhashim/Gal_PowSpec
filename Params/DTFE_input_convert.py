from __future__ import division
from numpy import *

x, y, z, lg10M, v_x, v_y, v_z  = loadtxt("/users/mhashim/lustre/mhashim/OxProj/Galaxy_Catalog/tmpcat_761_laigle.dat", unpack=True)
h = 0.704
f = open("/users/mhashim/lustre/mhashim/OxProj/Galaxy_Catalog/tmpcat_761_laigle_inputDTFE.dat", "w") 
f.writelines("%s\n" %(126361))
f.writelines("%f        %f      %f      %f      %f      %f\n" %(-50.0,50.0, -50.0,50.0, -50.0,50.0))
f.writelines("%f	%f	%f	%f\n" %(x[i]*h, y[i]*h, z[i]*h, 10.0**(lg10M[i] - 10.0)*h) for i in range(len(x)))
#f.writelines("%f       %f      %f      %f\n" %(x[i]*h, y[i]*h, z[i]*h, 10.0**(10.0 - 10.0)*h) for i in range(len(x)))
f.close()
