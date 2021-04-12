'''
Compute uncertainties of the effective phase noise and the theoretical expected growth rate.
Note: The uncertainties of the bunch length estimation are not included yet.
'''

import numpy as np
from math import *
import pickle as pkl
from my_functions import *

# Machine and beam parameters

betay = 73.81671646 # 73 m at CC2, 76 m at CC1
clight = 299792458 # m/s
Vcc = 1e6 # V
Eb = 270e9 # eV
gamma_0 = 287.8
beta_0 = np.sqrt(1 - 1/gamma_0**2)
circumference = 6911.5623 # m
frev = 299792458/circumference # Hz
numBunches = 4
gamma_0 = 287.7  # cmt it
f_RF = 400.789e6  # CC frequency in Hz

initial_bunch_length = False # if False the analysis is done using the average bunch length during each noise setting.

# Coast1-Setting1, Coast1-Setting2, Coast2-Setting1, Coast2-Setting2, Coast2-Setting2, Coast3-Setting1, Coast3-Setting2, Coast3-Setting3
PN_dBc = np.array([-122.75, -101.48, -115.22, -111.28, -111.03, -106.46, -101.48])
sigma_PN_dBc = np.array([0.47, 0.83, 0.61, 0.6, 0.83, 0.3, 0.76]) # uncertainty in the computation of the PN from the measured spectrum, averaged values over 1kHz
AN_dBc = np.array([-128.15, -115.21, -124.06, -115.71, -116.92, -112.73, -106.99])
sigma_AN_dBc = np.array([0.57, 0.66, 0.5, 0.4, 0.37, 0.3, 0.42 ])

# Load the computed effective phase noise for the noise levels mentioned above
if initial_bunch_length:
    effective_PN_dict = pkl.load(open('output/effective_PN_4bunches_MD5_initialBunchLengths.pkl.pkl', 'rb'))
else:
    effective_PN_dict = pkl.load(open('./output/effective_PN_4bunches_MD5_averagedBunchLengths.pkl', 'rb'))

# uncertainties of the measured PN and AN when converted in rad^2/Hz
# big differences on the uncertainties of the theoretical rates are introduced here
sigma_PN_rad = np.log(10)*ssb_2_dsb(PN_dBc)*sigma_PN_dBc/10
sigma_AN_rad = np.log(10)*ssb_2_dsb(AN_dBc)*sigma_AN_dBc/10


# uncertainty of the theoretical expected growth rate
# Load the initial bunch lengths for each coast settings
if initial_bunch_length:
    print('Use initial bunch lengths')
    with open('./bunchLengths_pickleFiles/initial_bunch_lengths_4sigmat_seconds.pkl', 'rb') as handle:
        b = pkl.load(handle)
else:
    print('Use averaged bunch lengths for each noise setting')
    with open('./bunchLengths_pickleFiles/averaged_bunch_lengths_4sigmat_seconds.pkl', 'rb') as handle:
        b = pkl.load(handle)
keys_b = list(b.keys())

sigma_phi_dict = {}
for key in keys_b: # iterate over the 4 bunches
    sigma_t = np.array(b[key])/4 # 1 sigma_t in seconds to be used in the function that converts to m
    sigma_z = bunch_length_time_to_m(sigma_t, clight)
    sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_RF)  # rms bunch length in rad
    sigma_phi_dict[key] = sigma_phi

sigma_theory_dey_dex_dict = {}
for key in keys_b:
    myC_PN = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], 'PN')
    myC_AN = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], 'AN')

    C_PN = 2*betay*(Vcc*frev/(2*Eb))**2*myC_PN
    C_AN = 4*betay * (Vcc * frev / (2 * Eb)) ** 2 * myC_AN

    sigma_dey_PN = C_PN*sigma_PN_rad
    sigma_dey_AN = C_AN * sigma_AN_rad



    sigma_theory_dey_dex_dict[key] = np.sqrt(sigma_dey_PN**2+sigma_dey_AN**2)

print('Uncertainty on the theoretical expected growths due to the uncertainty of the noise {}'.format(sigma_theory_dey_dex_dict[key]))

# uncertainty on the effective noise in rad^2/Hz
sigma_effective_PN_rad_dict = {}
for key in keys_b:
    myC_PN = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], 'PN')
    myC_AN = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], 'AN')

    C_PN = 2*betay * (Vcc * frev / (2 * Eb)) ** 2 * myC_PN
    C_AN = 4*betay * (Vcc * frev / (2 * Eb)) ** 2 * myC_AN

    sigma_dey_PN = C_PN * sigma_PN_rad
    sigma_dey_AN = C_AN * sigma_AN_rad

    sigma_effective_PN_rad_dict[key] = np.sqrt(sigma_dey_PN**2+sigma_dey_AN**2)/C_PN


# uncertainty on the effective noise in dbc/Hz
sigma_effective_PN_dBc_dict = {}

for key in keys_b:
    sigma_effective_PN_dBc_dict[key] = np.log10(e)*sigma_effective_PN_rad_dict[key]*10/ssb_2_dsb(np.array(effective_PN_dict[key])) # not sure for x10, maybe x5

print('Uncertainty for bunch 1 on the effective PN in dBc/Hz {}'.format(sigma_effective_PN_dBc_dict['bunch_1']))

save_bunch_lengths = True
if save_bunch_lengths:
    with open(f'./output/uncertainties_effective_PN_4bunches_MD5_initialBunchLength{str(initial_bunch_length)}.pkl', 'wb') as f: # 4sigma_t in s
        pkl.dump(sigma_effective_PN_dBc_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
    with open(f'./output/uncertainties_on_theoretical_rates_4bunches_MD5_initialBunchLength{str(initial_bunch_length)}.pkl', 'wb') as f: # 4sigma_t in s
        pkl.dump(sigma_theory_dey_dex_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

