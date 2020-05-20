import numpy as np
import numpy.ma as ma

def folding_reduced_chi_square(folded,foldedExpo,bins=10,weight=None,weightExpo=None):
    '''
    Calculate the chi square against uniform distribution for folded arrival times
    
    INPUT:
    folded: folded TOA in the dim of (phase,trial_period)
    foldedExpo: folded exposure time (phase, trial period)
    bins: number of phase bin
    weight: used to down-weight close pairs.
    weightExpo: used to down-weight exposure time (not used in this approach).
    
    OUTPUT:
    return chi square
    '''
    
    chi=np.zeros(len(folded))
    stdExposure=np.zeros(len(folded))

    step=1./bins
    if weight is None:
        weight=np.ones(len(folded[0]))
    if weightExpo is None:
        weightExpo=np.ones(len(foldedExpo[0]))
        
    for i in np.arange(len(folded)):
        #med=np.median(folded[i]) #place the median at the center of first bin
        hist,edge=np.histogram(folded[i],bins=np.arange(0,1+step,step),weights=weight)
        histExpo,edge=np.histogram(foldedExpo[i],bins=np.arange(0,1+step,step),normed=True,weights=weightExpo)
        
        exp=hist.sum()/histExpo.sum()
        histExpo=ma.array(histExpo,mask=histExpo<0.1)
        exphistExpo=histExpo*exp      
        
        chi[i]=\
        ((hist-exphistExpo)**2/exphistExpo).sum()/(histExpo.count()-1)#(bins-1)#
        
    return chi

