def stdgaminv(p,gam):
    import numpy as np
    import scipy.special as special
    abs_tol = 10**-3
    rel_tol = 10**-3

    x_max = 10**np.polyval([-0.009486738 ,0.03376901, 0.1151316, 0.2358172, 1.139717],log10(gam))
    max_iter = 200

    current_iter = 0
    while special.gammainc(x_max,gam)<max(p):
        current_iter+=1
        if current_iter>max_iter:
            print('Maximum specified probability is too high:{}'.format(max(p)))
        else:
            x_max *=1.5 
    x_min = 10**np.polyval([-0.0854665, 0.866249, -3.25511, 5.14328, -0.90924, -8.09135, 12.3393, -5.89628],log10(gam))
    current_iter = 0
    
    while special.gammainc(x_min,gam)>min(p):
        current_iter +=1

        if current_iter>max_iter:
            print('Minimum specified probability is too low:{}'.format(min(p)))
        else:
            x_max *=0.1

    n_check = 1000
    x_check = np.linspace(x_min,x_max,n_check)
    p_check = special.gammainc(x_check,gam)


    
    p_check, ind_u = np.unique(p_check,return_index = True)

    x_check