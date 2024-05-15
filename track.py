from wimps import get_drder, Fn2SI
from SRIM_derived_data import olivine
import numpy as np
# print(get_drder_one(10, 23,np.array([5]),np.array([1e-45]), frac=0,spin='SI'))

Er = np.logspace(-1,3,100)
A = olivine.atomic_masses
atomic_fraction = olivine.atomic_fractions
mx = 5 # GeV
sigma = 1e-45 #cm2
# # test Fn2SI
# Fn2 = Fn2SI(Er, A[0])
# np.save('test_results/Fn2.npy',Fn2)

# test drder 
drder = get_drder(Er, A, atomic_fraction,mx,sigma)
np.save('test_results/drder_5GeV_45.npy', drder)
print(drder.shape)