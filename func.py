import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import matplotlib.backends.backend_pdf as pdf
import scipy

def downweight_multiburst(seqt,threshTime=0.5):
    '''
    Given an input time serie
    INPUT:
        seqt: a sequence of TOA.
        
        OPTIONAL:
        threshTime: the thresh hold time for bursts to be considered correlated.Should be in same unit with seqt
        
    '''
    seqt=np.sort(seqt)
    dseqt=np.diff(seqt)
    weight=np.ones(len(seqt))
    cnt=1
    for i in np.arange(len(weight)):
        if i==(len(weight)-1):
            weight[i+1-cnt:i+1]/=cnt
            break

        elif dseqt[i]<threshTime:
            cnt+=1
        else:
            weight[i+1-cnt:i+1]/=cnt
            cnt=1
    return weight


def read_in_wiki_repeater(filename):
    gplist=np.loadtxt(filename,dtype=str)
    timesMult=[]
    names=[]
    times=[];name=''
    for j in np.arange(len(gplist)):
        a=gplist[j,1]
        if gplist[j,0]==name:
            times.append('20'+a[:2]+'-'+a[2:4]+'-'+a[4:6]+'T'+gplist[j,2])
        else:
            timesMult.append(times)
            times=[]
            names.append(name)
            times.append('20'+a[:2]+'-'+a[2:4]+'-'+a[4:6]+'T'+gplist[j,2])
            name=gplist[j,0]

    names.append(name)
    timesMult.append(times)
    timesMult=timesMult[1:]
    names=names[1:]
    return names,timesMult

def get_folded(seqt,foldPeriods):
    '''
    INPUT: 
        seqt: 1D TOA series
        foldPeriods: the 1D array of periods to fold
    OUTPUT:
        folded: the folded array with dim (periods, TOA)
    '''
    folded=((np.repeat(seqt[:,None],len(foldPeriods),axis=1)/foldPeriods)%1).T
    return folded


def h_test_repeaters(timesMult,names,foldPeriods,numBurstThresh=5,weightMethod=None)
    sid=86164.0905/3600/24

    htestOut=[]
    nameOut=[]
    foldedOut=[]
    sigmaOut=[]
    for i in np.arange(len(timesMult))[:]:
        if len(timesMult[i])>numBurstThresh:
            tTot=Time(timesMult[i])
            tTot=np.sort(tTot)
            seqt=np.array([(tTot[j]-tTot[0]).value for j in np.arange(len(tTot))])
            folded=get_folded(seqt,foldPeriods)

            htests=np.zeros((len(foldPeriods),3))
            if weightMethod is None:
                for j in np.arange(len(foldPeriods)):
                    htests[j]=h_test(folded[j],weights=None,max_harmonic=10)

            elif weightMethod=='day':
                weightData=downweight_multiburst(seqdt,threshTime=0.5)
                for j in np.arange(len(foldPeriods)):
                    htests[j]=h_test(folded[j],weights=weightData,max_harmonic=10)

            elif weightMethod=='cycle':
                foldTwice=np.concatenate((folded,folded+1),axis=-1)
                foldTwice=np.sort(foldTwice,axis=-1)
                gap=np.diff(foldTwice,axis=-1)
                offwindow=gap.max(-1)
                for j in np.arange(len(foldPeriods)):
                    threshTime=(1-offwindow[j])*foldPeriod[j]
                    weightData=downweight_multiburst(seqdt,threshTime=threshTime)
                    htests[j]=h_test(folded[j],weights=weightData,max_harmonic=10)

        trial=seqt[-1]/foldPeriods[0]
        chances=htests[:,-1]*trial;chances[chances>1]=1
        sigma=scipy.stats.norm.isf(chances/2.)

        htestOut.append(htests)
        nameOut.append(names[i]+'  %d bursts'%(len(seqt)))
        foldedOut.append(folded)
        sigmaOut.append(sigma)
    return nameOut,foldedOut,htestOut,sigmaOut

def plot_htest_out(nameOut,foldPeriods,htestOut,sigmaOut,figureName='repeater_htest.pdf'):
    page=pdf.PdfPages(figureName)
    for i in np.arange(len(htestOut)):
        htests=np.copy(htestOut[i])
        sigma=np.copy(sigmaOut[i])
        sel=np.logical_and(np.abs((foldPeriods/sid+0.5)%1-0.5)<0.1,foldPeriod/sid<10)

        fig,axes=plt.subplots(2,figsize=[10,8],sharex=False)
        axes[0].plot(foldPeriods[~sel],htests[:,0][~sel],'k.')
        axes[1].plot(foldPeriods[~sel],sigma[~sel],'k.')

        axes[0].set_xscale('log')
        axes[-1].set_xlabel('period day')
        axes[0].set_title(nameOut[i])
        axes[0].set_ylabel('h test')
        axes[1].set_ylabel('sigma')
        page.savefig(fig)
    page.close()

def multinomial_active_cycle_repeater(timesMult,names,foldPeriods,numBurstThresh=3,weightMethod=None)

    chancesOut=[]
    nameOut=[]
    foldedOut=[]
    numCycleOut=[]
    windowOut=[]
    sigmaOut=[]
    for i in np.arange(len(timesMult))[:]:
        if numEpoch[i]>numBurstThresh:

        #if len(timesMult[i])>10:

            tTot=Time(timesMult[i])
            tTot=np.sort(tTot)
            seqdt=np.array([(tTot[j]-tTot[0]).value for j in np.arange(len(tTot))])
            #foldPeriod=(np.arange(1./100,1./5,0.5/seqdt[-1]))**(-1)

            folded=((np.repeat(seqdt[:,None],len(foldPeriod),axis=1)/foldPeriod)%1).T
            foldTwice=np.concatenate((folded,folded+1),axis=-1)
            foldTwice=np.sort(foldTwice,axis=-1)
            gap=np.diff(foldTwice,axis=-1)
            offwindow=gap.max(-1)

            #count number of cycles
            dt=np.diff(seqdt)
            dtArr=np.repeat(dt[None,:],len(foldPeriod),axis=0)#dim periods,dt
            newCycle=(dtArr>=(offwindow*foldPeriod)[:,None])
            numCycles=np.count_nonzero(newCycle,axis=-1)+1
            chances=((1-offwindow)**(numCycles-1))
            chanceNew=np.copy(chances)*seqdt[-1]/foldPeriod[0] #times number of independent periods searched
            chanceNew[chanceNew>1]=1 #to avoid inf
            sigma=scipy.stats.norm.isf(chanceNew/2)


            nameOut.append(names[i])
            chancesOut.append(chances)
            foldedOut.append(foldTwice)
            numCycleOut.append(numCycles)
            windowOut.append(1-offwindow)
            sigmaOut.append(sigma)
            return nameOut,chancesOut,foldedOut,numCycleOut,windowOut,sigmaOut

def save_active_window_search_result(nameOut,chancesOut,foldedOut,numCycleOut,windowOut,sigmaOut,fileNameOut='period_candidates_active_window.txt'):
    with open(fileNameOut,'w') as f:
    f.write('#repeater, period candidate (d), number of cycles, active window fraction, \
    chances of coincidence (after times number of trials), sigma \n')
    for i in np.arange(len(nameOut)):
        sel=np.argsort(chancesOut[i])[0]
        f.write('{} {:.2f} {} {:.2f} {:.1e} {:.1f} \n'\
        .format(nameOut[i],foldPeriod[sel],numCycleOut[i][sel],windowOut[i][sel],chancesOut[i][sel],sigmaOut[i][sel]))


def plot_active_window_result(nameOut,foldPeriods,numCycleOut,windowOut,chancesOut,sigmaOut,figureName='repeater_activeWindow.pdf'):
    sid=86164.0905/3600/24
    sel=np.abs((foldPeriod/sid+0.5)%1-0.5)<0.01

    page=pdf.PdfPages(figureName)
    for i in np.arange(len(nameOut)):
        numCycles=numCycleOut[i]

        foldPeriodMask=np.copy(foldPeriod[~sel])
        pdg=np.copy(offwindow[~sel])


        fig,axes=plt.subplots(4,figsize=[10,8],sharex=True)
        axes[0].plot(foldPeriodMask,pdg,'k.')
        axes[0].set_xlabel('period day')
        axes[0].set_ylabel('off window percent')
        axes[0].set_title(nameOut[i])

        axes[1].plot(foldPeriodMask,numCycleOut[i][~sel],'k.')
        axes[1].set_ylabel('num of cycles')

        axes[2].plot(foldPeriodMask,chancesOut[i][~sel],'k.')
        axes[2].set_yscale('log')
        axes[2].set_ylabel('raw chances')

        axes[3].plot(foldPeriodMask,sigmaOut[i][~sel],'k.')
        axes[3].set_ylabel('sigma')
        axes[0].set_xscale('log')
        axes[-1].set_xlabel('period day')
        page.savefig(fig)
    page.close()
