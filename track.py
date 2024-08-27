from WIMPS.compute_wimps import get_drder, get_drder_one
from minerals_def import olivineObj
from compute_trackSpectra import get_wimps_Nbins, get_neutron_Nbins, get_wimps_dRdx_one,get_bkgNeutrino_dRdx_mol,get_neutron_dRdx_mol
from light_mediators import get_dRdE_nu
import numpy as np
from scipy.interpolate import interp1d

# binned rates

binwidth = 1
xt_range = np.logspace(-2,3,500)
mineral = olivineObj
mx = 5

sigma = 1e-45

# dRdx DM with SRIM

if False:
    xt_range = np.logspace(-2,3,500)
    mx = 5
    sigma  = np.array([1e-47])
    dRdx_dm = 0
    for i in range(4):
        A = olivineObj.atomic_masses[i]
        data_path =olivineObj.derived_data_path[i]
        fraction = olivineObj.atomic_fractions[i]

        dRdx_dm += get_wimps_dRdx_one(mx, sigma, xt_range, A, data_path)*fraction
    
    
    # np.save('/Users/szechingaudreyfung/PaleoBSMwithTRIM/test_results/interp_vals_olivine_SM_interp', dists)
    np.save('/Users/szechingaudreyfung/PaleoBSMwithTRIM/test_results/drdx_olivine_5GeV_47_interp', dRdx_dm)
    print(dRdx_dm.shape)

    

# neutron dRdx with srim
if True:
    xt_range = np.logspace(-2,3,500)
    # dRdx_neutron = 0
    # for i in range(4):
    #     data_path = olivineObj.derived_data_path[i]
    #     num_frac = olivineObj.number_fractions[i]

    #     dRdx_neutron += get_neutron_dRdx_one(xt_range, data_path, 153.31,1e-10)
    dRdx_neutron = get_neutron_dRdx_mol(xt_range, olivineObj,1e-10)
    
    np.save('/Users/szechingaudreyfung/PaleoBSMwithTRIM/test_results/drdx_olivine_neutron_benchmark_interp', dRdx_neutron)

# neutrino dRdx with srim
if True:
    xt_range = np.logspace(-2,3,500)
    mineral=olivineObj
    dRdx_bkgNu = get_bkgNeutrino_dRdx_mol(xt_range, mineral)
    np.save('/Users/szechingaudreyfung/PaleoBSMwithTRIM/test_results/drdx_olivine_bkgNu_interp', dRdx_bkgNu)


# test with Stopping power
if False:
    xt_range = np.logspace(-2,3,500)
    SP_files = ['OOli_SRIM.txt', 'SiOli_SRIM.txt', 'FeOli_SRIM.txt', 'MgOli_SRIM.txt']
    baum_fraction = np.array([0.41743753, 0.18319438, 0.14570679, 0.2536613 ])
    mx = 5
    sigma  = np.array([1e-47])
    dRdx_sp = 0
    for i in range(4):

        data_path =olivineObj.derived_data_path[i]
        fraction = olivineObj.atomic_fractions[i]
        A = olivineObj.atomic_masses[i]

        Erange = np.logspace(-3,3, 700)
        dRdE = get_drder_one(Erange, A,mx,sigma)

        dEdx = np.loadtxt(f'/Users/szechingaudreyfung/paleo/paleoSpec/Olivine/Ranges/{SP_files[i]}')
        dEdx = np.moveaxis(dEdx, 0,1)
        dEdx_interp = interp1d(dEdx[2]/10, dEdx[1]/baum_fraction[i]*10, fill_value=0, bounds_error=False)(xt_range)
        Erange_interp = interp1d(dEdx[2]/10, dEdx[0], fill_value=0, bounds_error=False)(xt_range)
        dRdE_interp = interp1d(Erange, dRdE, fill_value=0, bounds_error=False)(Erange_interp)

        dRdx_sp += dRdE_interp * dEdx_interp * fraction

    np.save('/Users/szechingaudreyfung/PaleoBSMwithTRIM/test_results/dRdx_5GeV_47_SP_olivine_benchmark.npy', dRdx_sp)
    print(dRdx_sp.shape)
