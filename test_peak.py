import numpy as np
from peakpressure.peakpressure import maxminest
from matplotlib import pyplot as plt

###########
# NOTE: seems to work with maxminest: read-in generated
given_series = np.loadtxt('test_data_gevrnd.dat', skiprows=0, usecols = (0,))
#given_series = np.loadtxt('test_data_normrnd.dat', skiprows=0, usecols = (0,))

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