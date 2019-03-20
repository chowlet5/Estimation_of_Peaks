def maxminest (record, dur_ratio = 1):

    import numpy as np
    import scipy.interpolate as interpolate
    import math
    #    if not record:


    n_cdf_pk =1000
    cdf_pk_min = 0.0005
    cdf_pk_max = 0.9995

    cdf_pk = np.linspace(cdf_pk_min,cdf_pk_max,n_cdf_pk)

    rsize = np.array(record).shape



    max_est = np.zeros((rsize[0],1))
    min_est = np.zeros((rsize[0],1))
    max_std = np.zeros((rsize[0],1))
    min_std = np.zeros((rsize[0],1))

    for i in np.arange(rsize[0]):
        x = record[i,:]
        n = len(x)

        mean_x = np.mean(x)
        std_x = np.std(x,ddof = 1)
        skew_x = np.sum(np.power(x-mean_x).3)/(n*std_x**3)
        X = x*np.sign(skew)

        sort_X = X.sort()

        mean_X = mean_x*np.sign(skew_x)
        std_X = std_x
        CDF_X = np.divide(np.arange(1,n+1),n+1)

        n_coarse = min([n,1000])

        CDF_coarse = np.linspace(1/(n_coarse+1),n_coarse/(n_coarse+1),n_coarse)

        f = interpolate.interp1d(CDF_X,sort_X)(CDF_coarse)
        X_coarse = f(CDF_coarse)

        mean_X_coarse = np.mean(X_coarse)
        std_X_coarse = np.std(X_coarse)

        gamma_min = 1
        gamma_max = 125
        n_gamma = 19
        n_start = 8

        gamma_list = np.logspace(math.log10(gamma_min),math.log10(gamma_max),n_gamma)
        gam_PPCC_list = zeros(gamma_list.shape)
        count = 0

        for j in np.arange(n_start,0,-1):
            count+=1

            s_gam_j = np.




     
