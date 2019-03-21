def stdnorminv(p):
    import numpy as np
    import scipy.special as special
    x = np.sqrt(2)*special.erfcinv(2*p)

    return x