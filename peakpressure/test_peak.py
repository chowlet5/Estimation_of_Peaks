import numpy as np
from peakpressure import maxminest
from matplotlib import pyplot as plt

###########
# NOTE: seems to work with maxminest: read-in generated
given_series = np.loadtxt('test_data_gevrnd.dat', skiprows=0, usecols = (0,))
#given_series = np.loadtxt('test_data_normrnd.dat', skiprows=0, usecols = (0,))

# ###########
# # NOTE: seems to work with maxminest: generate on the fly

# # start time
# start_time = 0.0
# # end time
# end_time = 10.0
# # steps 
# n_steps = 10000
# # time step
# delta_time = end_time / (n_steps-1)
# # time series
# # generate grid size vector (array) 1D
# time_series = np.arange(start_time, end_time + delta_time, delta_time)

# # frequency of the cosine
# cos_freq = 10
# # amplitude of the cosine
# cos_ampl = 1
# # series of the cosine
# cos_series = cos_ampl * np.cos(2*np.pi * cos_freq * time_series)

# # amplitude of the constant
# const_ampl = 10
# # series of the constant
# const_series = const_ampl * np.ones(len(time_series))

# # random signal 
# # assuming nomarl distribution
# # with given mean m = 0 and standard deviation std = 0.25
# rand_m = 0.0
# rand_std = 0.25
# # series of the random
# # NOTE seems to work properly also in Python as in Matlab after updating original script
# # rand_series = np.random.normal(rand_m, rand_std, len(time_series))
# k = 0.0
# sigma = rand_std
# mu = rand_std
# # series of gev -> gumbel
# rand_series = np.random.gumbel(rand_m, rand_std, len(time_series))

# # coefs -> weighting factors for the respective series of signals
# coef_signal1 = 1
# coef_signal2 = 0.25
# coef_signal3 = 1
# superposed_series = coef_signal1 * const_series + coef_signal2 * cos_series + coef_signal3 * rand_series

# ###########
# # NOTE: best to store randomly generated signal otherwise test values will alter during each new generation
# given_series = superposed_series

###########
# print results
result = maxminest(given_series)
print(result)

'''
RESULTS FOR 
    test_data_gevrnd.dat:
2.36535988
-0.71199031 
0.06569003
0.00299395

    test_data_normrnd.dat:
1.13118079
-1.0013704
0.00886154
0.00542217
'''
maxv = result[0][0][0]
minv = result[1][0][0]

###########

plt.figure()
x_series = np.arange(0.0, len(given_series), 1.0)
plt.plot(x_series, given_series)
plt.ylabel('Amplitude')
plt.xlabel('Time [s]')
plt.hlines([maxv, minv], x_series[0], x_series[-1])
plt.grid(True)
plt.show()