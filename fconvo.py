import numpy as np

#***************************************************************
#***************************************************************
#***************************************************************

def rebin(wv_in, fl_in, wv_out):
    """
    Rebins a spectrum onto a new array using linear interpolation

    wv_in and fl_in are the wavelength and flux of the original spectrum;
    wv_new is the wavelength to which to rebin

    returns a tuple of two vectors: 
    fl_new is the rebinned flux and 
    flag is the vector of flags: True if at the final location in wv_new the flux
    is actually defined, and False if not because the final wavelength baseline is 
    larger than the original

    linear interpolation is used
    this function assumes that both wavelength arrays are sorted and returns zero
    otherwise
    this function does not check that the wavelength arrays are in the same units,
    this is the user's responsibility
    """
    if (np.sum((np.sort(wv_in)-wv_in)**2)>1e-5): 
        print("input array not sorted")
        return(0)
    if (np.sum((np.sort(wv_out)-wv_out)**2)>1e-5): 
        print("output array not sorted")
        return(0)
    flag=np.ones(np.size(wv_out))
    flag[(wv_out<wv_in[0])]=0.
    flag[(wv_out>wv_in[np.size(wv_in)-1])]=0.

    # the middle part where the interpolation actually happens
    wv2=wv_out[(flag==1)]
    # and the other parts
    wv1=wv_out[(wv_out<wv_in[0])]
    wv3=wv_out[(wv_out>wv_in[np.size(wv_in)-1])]
    # for every new wavelength element, we need to identify the index of the
    # old wavelength element preceding this one and one that's exceeding this one
    tnum=np.zeros(np.size(wv2))
    uu=0
    tt=0
    while ((tt<np.sum(flag)) & (uu<np.size(wv_in)-1)):
        if ((wv2[tt]>=wv_in[uu]) & (wv2[tt]<=wv_in[uu+1])):
            tnum[tt]=uu
            tt=tt+1
            #print("I am here 1",tt,uu)
        else:
            uu=uu+1    
            #print("I am here 2",tt,uu)
    tnum=tnum.astype(int)
    fl_new2=fl_in[tnum]+(fl_in[tnum+1]-fl_in[tnum])/(wv_in[tnum+1]-wv_in[tnum])*(wv2-wv_in[tnum])
    fl_new=np.concatenate((0.*wv1,fl_new2,0.*wv3))
    return(fl_new,flag)

#***************************************************************
#***************************************************************
#***************************************************************

def fnu_ave(wv, fnu, wv_tran, tran, tflag=0, verbose=1):
    """
    This function averages a spectrum using a filter curve using equations from
    https://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/

    variables wv and fnu are wavelength and f_nu of the input spectrum, sorted by ascending wavelength
    wv_tran and tran are wavelength and response of the transmission curve, sorted by ascending wavelength
    the function assumes that wv and wv_tran are in the same units
    
    last one is the transmission flag:
    tflag=1 for photon counter or a quantum efficiency transmission curve
    tflag=0 for energy counter or bolometer
    
    returns filter-averaged fnu in the same units as input

    For SDSS, use tflag=1: http://classic.sdss.org/dr2/products/general/edr_html/node25.html 
    (the svo2 filter collection appears to give the wrong information about the SDSS filters)
    
    For WISE, Spitzer and Herschel filters use tflag=0. 
    WISE: http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html
    Spitzer IRAC: https://ui.adsabs.harvard.edu/abs/2005PASP..117..978R/abstract equation 9
    Spitzer MIPS: http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
    Herschel SPIRE: http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Herschel/SPIRE.PSW&&mode=browse&gname=Herschel&gname2=SPIRE#filter
    
    For HST WFC3 filters use tflag=0
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=HST/WFC3_IR.F105W&&mode=browse&gname=HST&gname2=WFC3_IR#filter
    For PAN-STARRS, use tflag=0
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=PAN-STARRS/PS1.g&&mode=browse&gname=PAN-STARRS&gname2=PS1#filter
    
    """
    if ((tflag!=0) & (tflag!=1)):
        print('unknown type of transmission curve')
        return(0)
    # The primary method is to rebin the spectrum onto the filter curve
    #from scipy.interpolate import interp1d
    # rebin the spectrum onto the filter curve: 
    #interp_lin = interp1d(wv, fnu)
    #_temp=interp_lin(wv_tran)
    (_temp,_tempflag)=rebin(wv,fnu,wv_tran)
    _numer=_temp*wv_tran**(tflag-2.)*tran
    _denom=wv_tran**(tflag-2.)*tran
    _ind=np.arange(0,len(_numer)-1)
    fnuave=np.sum((_numer[_ind+1]+_numer[_ind])*(wv_tran[_ind+1]-wv_tran[_ind]))/np.sum((_denom[_ind+1]
                                +_denom[_ind])*(wv_tran[_ind+1]-wv_tran[_ind]))
    
    if (verbose):
        # OK, this doesn't work if one of the curves is beyond the limit coverage of the other one... 
        #interp_lin = interp1d(wv_tran, tran)
        #_temp=interp_lin(wv) 
        (_temp,_tempflag)=rebin(wv_tran,tran,wv)
        # beyond the filter curve data, set it to 0 just in case
        _ind=((wv<wv_tran[0]) & (wv>wv_tran[len(wv_tran)-1]))
        _temp[_ind]=0.
        _numer=fnu*wv**(tflag-2.)*_temp
        _denom=wv**(tflag-2.)*_temp
        _ind=np.arange(0,len(_numer)-1)
        other_fnuave=np.sum((_numer[_ind+1]+_numer[_ind])*(wv[_ind+1]-wv[_ind]))/np.sum((_denom[_ind+1]
                                    +_denom[_ind])*(wv[_ind+1]-wv[_ind]))
        # the relative difference between the two results
        print('relative difference in fnu calculated by diff methods=', np.abs(other_fnuave-fnuave)/fnuave)
        
    return(fnuave)

#***************************************************************
#***************************************************************
#***************************************************************

def flam_ave(wv, flam, wv_tran, tran, tflag=0, verbose=1):
    """
    This function averages a spectrum using a filter curve using equations from
    https://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/

    variables wv and flam are wavelength and f_lambda of the input spectrum, sorted by ascending wavelength
    wv_tran and tran are wavelength and response of the transmission curve, sorted by ascending wavelength
    the function assumes that wv and wv_tran are in the same units
    
    last one is the transmission flag:
    tflag=1 for photon counter or a quantum efficiency transmission curve
    tflag=0 for energy counter or bolometer
    
    returns filter-averaged flam in the same units as input

    For SDSS, use tflag=1: http://classic.sdss.org/dr2/products/general/edr_html/node25.html 
    (the svo2 filter collection appears to give the wrong information about the SDSS filters)
    
    For WISE, Spitzer and Herschel filters use tflag=0. 
    WISE: http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html
    Spitzer IRAC: https://ui.adsabs.harvard.edu/abs/2005PASP..117..978R/abstract equation 9
    Spitzer MIPS: http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php
    Herschel SPIRE: http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Herschel/SPIRE.PSW&&mode=browse&gname=Herschel&gname2=SPIRE#filter
    
    For HST WFC3 filters use tflag=0
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=HST/WFC3_IR.F105W&&mode=browse&gname=HST&gname2=WFC3_IR#filter
    For PAN-STARRS, use tflag=0
    http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=PAN-STARRS/PS1.g&&mode=browse&gname=PAN-STARRS&gname2=PS1#filter
    
    """
    if ((tflag!=0) & (tflag!=1)):
        print('unknown type of transmission curve')
        return(0)
    # The primary method is to rebin the spectrum onto the filter curve
    #from scipy.interpolate import interp1d
    # rebin the spectrum onto the filter curve: 
    #interp_lin = interp1d(wv, flam)
    #_temp=interp_lin(wv_tran)
    (_temp,_tempflag)=rebin(wv,flam,wv_tran)
    _numer=_temp*wv_tran**tflag*tran
    _denom=wv_tran**tflag*tran
    _ind=np.arange(0,len(_numer)-1)
    flamave=np.sum((_numer[_ind+1]+_numer[_ind])*(wv_tran[_ind+1]-wv_tran[_ind]))/np.sum((_denom[_ind+1]
                                +_denom[_ind])*(wv_tran[_ind+1]-wv_tran[_ind]))
    
    if (verbose):
        # OK, this doesn't work if one of the curves is beyond the limit coverage of the other one... 
        #interp_lin = interp1d(wv_tran, tran)
        #_temp=interp_lin(wv) 
        (_temp,_tempflag)=rebin(wv_tran,tran,wv)
        # beyond the filter curve data, set it to 0 just in case
        _ind=((wv<wv_tran[0]) & (wv>wv_tran[len(wv_tran)-1]))
        _temp[_ind]=0.
        _numer=flam*wv**tflag*_temp
        _denom=wv**tflag*_temp
        _ind=np.arange(0,len(_numer)-1)
        other_flamave=np.sum((_numer[_ind+1]+_numer[_ind])*(wv[_ind+1]-wv[_ind]))/np.sum((_denom[_ind+1]
                                    +_denom[_ind])*(wv[_ind+1]-wv[_ind]))
        # the relative difference between the two results
        print('relative difference in flam calculated by diff methods=', np.abs(other_flamave-flamave)/flamave)
        
    return(flamave)

