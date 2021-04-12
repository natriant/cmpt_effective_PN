import numpy as np
import pickle as pkl
from my_functions import *

# Machine and beam parameters
betay = 73  # m in CC2, ~76 in CC1 (MAD-X)
Vcc = 1e6  # V
frev = 43.45e3  # Hz
Eb = 270e9  # eV
sigma_z = 0.138  # rms bunch length, m
beta_0 = 0.999999  # cmpt it from the rest
gamma_0 = 287.7  # cmt it
clight = 299792458  # light speed in meters/second
f_RF = 400.789e6  # CC frequency in Hz

# compute correction factor due to bunch length. Different fro each bunch according to the MD conditions

sigma_phi = bunch_length_m_to_rad(sigma_z, clight, f_RF)  # rms bunch length in rad

# Compute the expected growth rates from the measured AN and PN independently
PSD_PN_list = [-122.75, -101.48, -115.22, -111.28, -111.03, -106.46, -101.48]  # dBc/Hz
PSD_AN_list = [-128.15, -115.21, -124.06, -115.71, -116.92, -112.73, -106.99]

PSD_PN_list_radHz = [ssb_2_dsb(x) for x in PSD_PN_list]  # rad^2/Hz
PSD_AN_list_radHz = [ssb_2_dsb(x) for x in PSD_AN_list]  # rad^2/Hz

# PSD_phi_list = [-117]
# PSD_list = [ssb_2_dsb(x) for x in PSD_phi_list] # rad^2/Hz
print('PSD initial {} rad^2/Hz'.format(PSD_PN_list))
print('PSD initial {} rad^2/Hz'.format(PSD_AN_list))

dey_AN_list = []  # m/s geometric emittance
dey_PN_list = []

for PSD in PSD_PN_list_radHz:
    my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='PN')
    print('correction factor due to bunch length = {}'.format(my_C))
    dey_PN_list.append(emit_growth_phase_noise(betay, Vcc, frev, Eb, my_C, PSD, True))  # true for one-sided PSD

for PSD in PSD_AN_list_radHz:
    my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='AN')
    print('correction factor due to bunch length = {}'.format(my_C))
    dey_AN_list.append(emit_growth_amplitude_noise(betay, Vcc, frev, Eb, my_C, PSD, True))  # true for one-sided PSD

# add the growth rates as AN and PN can be considerate independent
dey_total_list = []
for index, k in enumerate(dey_AN_list):
    dey_total_list.append(dey_AN_list[index] + dey_PN_list[index])

print('dey/dt geometric in nm/s: {}'.format(np.array(dey_total_list) * 1e9))
#print(('dey/dt normalised in nm/s: {}'.format(np.array(dey_total_list) * 1e9 * beta_0 * gamma_0)))
#print(('dey/dt normalised in um/h: {}'.format(np.array(dey_total_list) * 1e9 * beta_0 * gamma_0 * 1e-3 * 3600)))



# For one value of bunch length
# Compute the required phase noise (effective noise) to obtain these rates
my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='PN')
psd_list_2_init = cmpt_phase_noise_from_growth_rate(betay, Vcc, frev, Eb, my_C, np.array(dey_total_list), True)
psd_list_2 = [dsb_2_ssb(x) for x in psd_list_2_init]
print(psd_list_2)

# sanity check
dey_PN_list_sanity_test = []
for PSD in psd_list_2_init:
    my_C = cmpt_bunch_length_correction_factor(sigma_phi, noise_type='PN')
    dey_PN_list_sanity_test.append(emit_growth_phase_noise(betay, Vcc, frev, Eb, my_C, PSD, True))  # true for one-sided PSD

print('sanity check {}'.format(np.array(dey_PN_list_sanity_test)*1e9))
