import numpy as np
def isprime(n):
    if not (isinteger(n) and n > 0):
        msg = "Argument must be a natural number"
        raise ValueError(msg)
    if n == 1:
        return False
    elif n == 2:
        return True
    else:
        term1 = np.floor((np.math.factorial(n-1)+1)/n)
        term2 = ((np.math.factorial(n-1)-n+1)/n)
        g_n = (n-2)*(term1-term2)
        if g_n == n:
            return True
        else:
            return False
          
