A trial pipeline to look for long time scale periodicity on CHIME data. With three searching functions as described below.

# minimize active window and apply multinomial test. Specifically for repeater with small number of bursts.


# chi square test in binned folded profile: histChisquare.py
Not optimal in search, but intuitive, easy to add extensions to correct for non-poisson effect and easy to verify with multinomial tests when the sample size is small.

# h test:  htest.py 
apply searches in fourier space, more robust (would usually get a higher S/N then the chi square test above, but less intuitive to add extensions.  
