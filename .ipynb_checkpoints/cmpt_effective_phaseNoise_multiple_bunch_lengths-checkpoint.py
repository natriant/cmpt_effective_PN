import numpy as np
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

# compute correction factor due to bunch length. Different fro each bunch according to the MD conditions
# load initial bunch lengths for each setting of the MD, the order is:
# [coast1-setting1, coast1-setting2, coast2-setting1, coast2-setting2, coast3-setting1, coast3-setting2, coast3-setting3]
# load from  '/eos/user/n/natriant/2020/Longitudinal_plane_analysis_MD5_2018/ABWLM/'
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

#sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_RF)  # rms bunch length in rad

# Compute the expected growth rates from the measured AN and PN independently
PSD_PN_list = [-122.75, -101.48, -115.22, -111.28, -111.03, -106.46, -101.48]  # dBc/Hz
PSD_AN_list = [-128.15, -115.21, -124.06, -115.71, -116.92, -112.73, -106.99]

PSD_PN_list_radHz = [ssb_2_dsb(x) for x in PSD_PN_list]  # rad^2/Hz
PSD_AN_list_radHz = [ssb_2_dsb(x) for x in PSD_AN_list]  # rad^2/Hz

# PSD_phi_list = [-117]
# PSD_list = [ssb_2_dsb(x) for x in PSD_phi_list] # rad^2/Hz
print('PSD initial {} rad^2/Hz'.format(PSD_PN_list))
print('PSD initial {} rad^2/Hz'.format(PSD_AN_list))

dey_AN_list = {}  # m/s geometric emittance
dey_PN_list = {}

for key in keys_b: # iterate over the bunches
    dey_AN_list[key] = []
    dey_PN_list[key] = []
    # for PSD in PSD_PN_list_radHz:
    my_C_pn = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], noise_type='PN')  # type: array
    #print('correction factor due to bunch length = {}'.format(my_C_pn))
    dey_PN_list[key].append(emit_growth_phase_noise(betay, Vcc, frev, Eb, my_C_pn, np.array(PSD_PN_list_radHz), True))  # true for one-sided PSD

    # for PSD in PSD_AN_list_radHz:
    my_C_an = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='AN')  # type: array
    #print('correction factor due to bunch length = {}'.format(my_C_an))
    dey_AN_list[key].append(emit_growth_amplitude_noise(betay, Vcc, frev, Eb, my_C_an, np.array(PSD_AN_list_radHz), True))  # true for one-sided PSD

# add the growth rates as AN and PN can be considerate independent
dey_total_dict = {}
for key in keys_b:
    dey_total_dict[key] = np.array(dey_PN_list[key][0]+dey_AN_list[key][0])

print('dey/dt geometric in nm/s: {}'.format(np.array(dey_total_dict['bunch_1']) * 1e9))
#print(('dey/dt normalised in nm/s: {}'.format(np.array(dey_total_list) * 1e9 * beta_0 * gamma_0)))
#print(('dey/dt normalised in um/h: {}'.format(np.array(dey_total_list) * 1e9 * beta_0 * gamma_0 * 1e-3 * 3600)))


'''
# For one value of bunch length
# Compute the required phase noise (effective noise) to obtain these rates
my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='PN')
psd_list_2_init = cmpt_phase_noise_from_growth_rate(betay, Vcc, frev, Eb, my_C, np.array(dey_total_list), True)
psd_list_2 = [dsb_2_ssb(x) for x in psd_list_2_init]
print(psd_list_2)
'''
# For multiple bunch length values
psd_list_2_dict = {}
for key in keys_b:  # iterate over bunches. Each array contains different settings
    my_C = cmpt_bunch_length_correction_factor(sigma_phi_dict[key], noise_type='PN')
    psd_list_2_init = cmpt_phase_noise_from_growth_rate(betay, Vcc, frev, Eb, my_C, np.array(dey_total_dict[key]), True)
    psd_list_2_dict[key] = [dsb_2_ssb(x) for x in psd_list_2_init]

save_bunch_lengths = False

if save_bunch_lengths:
    with open(f'./output/effective_PN_4bunches_MD5_initialBunchLengths{str(initial_bunch_length)}.pkl', 'wb') as f: # 4sigma_t in s
        pkl.dump(psd_list_2_dict, f, protocol=pkl.HIGHEST_PROTOCOL)


quit()
# sanity check
dey_PN_list_sanity_test = []
for PSD in psd_list_2_init:
    my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='PN')
    dey_PN_list_sanity_test.append(emit_growth_phase_noise(betay, Vcc, frev, Eb, my_C, PSD, True))  # true for one-sided PSD

print('sanity check {}'.format(np.array(dey_PN_list_sanity_test)*1e9))
