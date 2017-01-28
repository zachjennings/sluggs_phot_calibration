"""
Useful functions for analyzing optimal apertures and calculating aperture corrections in data.

All functions assume data array is in standard format:
[RA,DEC,(mag/merr pairs for all apertures),fwhm in pix, ellipticity, photometry flags]
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def make_diff_arr(data,n_apers,cor_aper=None):
    """
    Create an array storing the differences between measurements in every aperture with
    a specified "corrected" apertures.
    
    data: 1d numpy array containing standard photometry data format
    n_apers: int, number of measured apertures in our photometry
    cor_aper: index of aperture that we wish to compare our magnitudes to
    """
    n_paers = data.shape[0]
    if cor_aper is None:
        cor_aper = n_apers - 1
        
    mags = data[:,2:n_apers+2]
    diff_arr = np.zeros((data.shape[0],n_apers))
    for i in range(n_apers):
        diff_arr[:,i] = mags[:,cor_aper] - mags[:,i]
        
    return diff_arr

def make_merr_arr(data,n_apers):
    return data[:,n_apers+2:(2.*n_apers)+2]

def calc_aper_cor(data,fwhm_max=4.5,fwhm_min=4.,n_apers=18,meas_aper=3,\
                  cor_aper=13,sat_cut=17.,faint_cut=20.,cir_cut=1.,flag_cut=1):
    """
    Function to calculate aperture corrections.       
                  
    data: 1d numpy array, standard format
    
    photometric cut parameters:
        fwhm_max: (float), maximum fwhm in pix to be considered a point source
        fwhm_min: (float), min pix to be considered a point source
    """
    n_apers = n_apers
    ra = data[:,0]
    dec = data[:,1]
    mags = data[:,2:n_apers+2]
    merrs = data[:,n_apers+2:(2.*n_apers)+2]
    flags = data[:,-1]
    fwhm = data[:,-3]
    ellip = data[:,-2]
    
    good = (flags < flag_cut) & (ellip < cir_cut) & (mags[:,meas_aper] > sat_cut) \
        & (mags[:,meas_aper] < faint_cut) & \
         (fwhm < fwhm_max)  & (fwhm > fwhm_min)
            
    aper_cor = np.median(mags[good,meas_aper] - mags[good,cor_aper])
    
    return (aper_cor,good)

def make_merr_plot(data,apers,figsize=(10,5),save=False,fname='merr_plot.png',dpi=300,**args):
    merrs = make_merr_arr(data,apers.shape[0])
    fig,ax = plt.subplots(figsize=figsize)
    ax.plot(apers,np.median(merrs,axis=0),**args)
    
    if save:
        fig.savefig(fname,dpi=dpi)
        
def make_cog_plot(data,apers,figsize=(10,5),cor_aper=10,save=False,fname='merr_plot.png',dpi=300,ylim=(-1.0,0.2),meas_aper=1):
    diff_arr = make_diff_arr(data,apers.shape[0],cor_aper=cor_aper)
    fig,ax = plt.subplots(figsize=figsize)
    medians = np.median(diff_arr,axis=0)
    
    ax.plot(apers,medians)
    ax.set_ylim(ylim)
    ax.axvline(x=apers[meas_aper],color='r')
    ax.axvline(x=apers[cor_aper],color='b')
    
    
    if save:
        fig.savefig(fname,dpi=dpi)
        
def make_diff_plot(data,good,meas_aper=2,cor_aper=10,min_plot=13,max_plot=25,n_apers=14,xlim=(14,25),ylim=(-2,2)):
    diff_arr = make_diff_arr(data,n_apers=n_apers,cor_aper=cor_aper)
    fig,ax = plt.subplots(figsize=(10,5))
    ax.scatter(data[:,meas_aper+2],diff_arr[:,meas_aper],s=0.1)
    ax.scatter(data[good,meas_aper+2],diff_arr[good,meas_aper],s=10,color='b')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
def write_reg(ra,dec,size='2"',color='green',shape='circle',filename='regions.reg'):
    f = open(filename,'w')
    f.write('fk5\n')
    for i,j in zip(ra,dec):
        f.write(shape+' '+str(i)+' '+str(j)+' '+size+' # color='+color+'\n')
    f.close()     