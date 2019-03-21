def stdnormcdf(x):
    import scipy.special as special
    import math
    p = 0.5 *special.erfc(-x/math.sqrt(2))
    return p