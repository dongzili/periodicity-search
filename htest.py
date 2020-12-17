import numpy.ma as ma
def h_test(phases, max_harmonic=20):
    """Apply the H test for uniformity on [0,1).
    The H test is an extension of the Z_m^2 or Rayleigh tests for
    uniformity on the circle. These tests estimate the Fourier coefficients
    of the distribution and compare them with the values predicted for
    a uniform distribution, but they require the user to specify the number
    of harmonics to use. The H test automatically selects the number of
    harmonics to use based on the data. The returned statistic, H, has mean
    and standard deviation approximately 2.51, but its significance should
    be evaluated with the routine h_fpp. This is done automatically in this
    routine.
    Arguments
    ---------
    events : array-like
        events should consist of an array of values to be interpreted as
        values modulo 1. These events will be tested for statistically
        significant deviations from uniformity.
    Returns
    -------
    H : float
        The raw score. Larger numbers indicate more non-uniformity.
    M : int
        The number of harmonics that give the most significant deviation
        from uniformity.
    fpp : float
        The probability of an H score this large arising from sampling a
        uniform distribution.
    Reference
    ---------
    de Jager, O. C., Swanepoel, J. W. H, and Raubenheimer, B. C., "A
    powerful test for weak periodic signals of unknown light curve shape
    in sparse data", Astron. Astrophys. 221, 180-190, 1989.
    Updated false alarm rate  to match Jager, Busching 2010
    ---------
    credit to Robert Archibald
    """
    if len(phases)==0:
        H=0
        M=0
        fpp=1
    else:
        ev = np.reshape(phases, (-1,))
        cs = np.sum(np.exp(2.j*np.pi*np.arange(1,max_harmonic+1)*ev[:,None]),axis=0)/len(ev)
        Zm2 = 2*len(ev)*np.cumsum(np.abs(cs)**2)
        Hcand = (Zm2 - 4*np.arange(1,max_harmonic+1) + 4)
        M = np.argmax(Hcand)+1
        H = Hcand[M-1]
        fpp =np.exp(-0.4*H)
    return (H, M, fpp)

def h_test_weight(phases, max_harmonic=20,weights=None):
    """Apply the H test for uniformity on [0,1). with weight of each data points
    The H test is an extension of the Z_m^2 or Rayleigh tests for
    uniformity on the circle. These tests estimate the Fourier coefficients
    of the distribution and compare them with the values predicted for
    a uniform distribution, but they require the user to specify the number
    of harmonics to use. The H test automatically selects the number of
    harmonics to use based on the data. The returned statistic, H, has mean
    and standard deviation approximately 2.51, but its significance should
    be evaluated with the routine h_fpp. This is done automatically in this
    routine.
    ---------
    version: 2020, support weights for the data points. Edited by Dongzi
    Original script: Robert Archbald
    ----------    
    Arguments
    ---------
    phases : array-like
        phases of events should consist of an array of values to be interpreted as
        values modulo 1. These events will be tested for statistically
        significant deviations from uniformity.
    max_harmonic=20: the maximum harmonic to be searched
    weights: default None. The weight of each event.
    
    Returns
    -------
    H : float
        The raw score. Larger numbers indicate more non-uniformity.
    M : int
        The number of harmonics that give the most significant deviation
        from uniformity.
    fpp : float
        The probability of an H score this large arising from sampling a
        uniform distribution.
    Reference
    ---------
    de Jager, O. C., Swanepoel, J. W. H, and Raubenheimer, B. C., "A
    powerful test for weak periodic signals of unknown light curve shape
    in sparse data", Astron. Astrophys. 221, 180-190, 1989.
    Updated false alarm rate  to match Jager, Busching 2010
    """
    if weights is None:
        weights=np.ones(len(phases))
        
    if len(phases)==0:
        H=0
        M=0
        fpp=1
    else:
        ev = np.reshape(phases, (-1,))
        cs = np.sum(np.exp(2.j*np.pi*np.arange(1,max_harmonic+1)*ev[:,None])*weights[:,None],axis=0)/weights.sum()
        Zm2 = 2*weights.sum()*np.cumsum(np.abs(cs)**2)
        Hcand = (Zm2 - 4*np.arange(1,max_harmonic+1) + 4)
        M = np.argmax(Hcand)+1 #optimal harmonic number
        H = Hcand[M-1] #raw score
        fpp =np.exp(-0.4*H) #probability from uniform
    return (H, M, fpp)
