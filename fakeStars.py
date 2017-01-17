import numpy as np
import subprocess as sp
import os
import pickle
import astropy.coordinates as coords
import astropy.units as u
#from sklearn.svm import SVC as svc
#from sklearn.preprocessing import StandardScaler
#from sklearn.cross_validation import StratifiedShuffleSplit
#from sklearn.grid_search import GridSearchCV
from astropy.wcs import WCS
from astropy.io import fits
import seWrapper as se

from pyraf.iraf import artdata


class fakeStarTests(object):
    def __init__(self,image,phot_catalog,pix_scale=0.2,gain=1.0,rd_noise=0.0,\
        zpt=30.,seeing=3.5,config_file='foo.sex',write_true_phot=False,exptime=1.,background=0.):
        #self.runner = seRunner()
        self.phot = phot_catalog
        self.name = image
        hdu = fits.open(image)
        self.header = hdu[0].header
        self.pix_scale=pix_scale
        self.gain=gain
        self.rd_noise=rd_noise
        self.zpt = zpt
        self.seeing = seeing
        self.config_file=config_file
        self.wcs = WCS(image)
        self.write_true_phot = write_true_phot
        self.exptime = exptime
        self.background = background

        if self.write_true_phot:
            self.true_mags,self.true_x,self.true_y = self.unpackPhotCatalog(self.phot)

    def unpackPhotCatalog(self,cat):
        true_x,true_y = self.wcs.wcs_world2pix(cat.coords.ra.deg,cat.coords.dec.deg,1)
        true_mags = cat.mags[cat.mags.keys()[0]]
        return true_mags,true_x,true_y

    def placeStars(self,name='',output='',objects=''):
        '''
        Use pyraf artdata.mkobjects to place the artificial stars in the image
        '''
        artdata.mkobjects(input=name,output=output,objects=objects,\
                        magzero=self.zpt,gain=self.gain,rdnoise=self.rd_noise,radius=self.seeing,background=self.background,exptime=self.exptime)


    def makeStarList(self,name,n_stars=1000,min_mag=24,max_mag=28,x_range=(0,1000),y_range=(0,1000),write_true_phot=False):
        '''
        Use pyraf artdata.starlist to make an input starlist
        '''
        #artdata.starlist(name,nstars=n_stars,xmax=self.header['naxis1'],ymax=self.header['naxis2'],\
        #    minmag=min_mag,maxmag=max_mag,luminosity='uniform',lseed='INDEF')

        mags = np.random.uniform(low=min_mag,high=max_mag,size=n_stars)
        x = np.random.uniform(low=x_range[0],high=x_range[1],size=n_stars)
        y = np.random.uniform(low=y_range[0],high=y_range[1],size=n_stars)

        fake_cat = np.vstack([x,y,mags])

        #If we're also writing out the true clusters in this step, need to re-make
        #the starlist so that it includes the photometry.
        if self.write_true_phot:
            full_mags = np.concatenate([self.true_mags,fake_cat[:,2]])
            full_x = np.concatenate([self.true_x,fake_cat[:,0]])
            full_y = np.concatenate([self.true_y,fake_cat[:,1]])
        
        else:
            full_mags = mags
            full_x = x
            full_y = y


        write_array = np.vstack([full_x,full_y,full_mags]).T

        sp.call('rm '+name,shell=True)
        with open(name,"w+") as the_file:
            for i in write_array:
                line = np.array2string(i,max_line_width=np.inf).replace('[','').replace(']','').strip()

                the_file.write(line + '\n')

        #return the fake catalogs so that we know what to check against
        return (x,y,mags)

    def genFakeCatalogs(self,n_iter=1,min_mag=22.,max_mag=26.,n_stars=1000,overwrite=False,xmax=None,ymax=None,xmin=0.,ymin=0.,match_rad=0.5*u.arcsec):
        '''
        Perform photometry on the fake catalogs
        '''
        if xmax is None:
            xmax = self.header['naxis1']
            
        if ymax is None:
            ymax = self.header['naxis2']
             
        
        if overwrite or not hasattr(self,'fake_cats'):
            self.fake_cats = dict()
            self.fake_cats_x = dict()
            self.fake_cats_y = dict()
            ini = True
        else:
            ini = False

        for i in range(n_iter):
            starlist_name = self.name+'.starlist'+str(i)
            fake_image = self.name+'.fake'+str(i)+'.fits'
            recovered_cat_name = fake_image+'.cat'
            sp.call('rm '+starlist_name,shell=True)
            sp.call('rm '+fake_image,shell=True)


            starlist_x,starlist_y,starlist_mags = self.makeStarList(starlist_name,min_mag=min_mag,max_mag=max_mag,n_stars=n_stars,\
                x_range=(xmin,xmax),y_range=(ymin,ymax))
            starlist_mags={'input':starlist_mags}
            starlist_fwhm = {'input':np.zeros(starlist_x.size)+99.999}
            starlist_ellipticity = {'input':np.zeros(starlist_x.size)+99.999}

            self.placeStars(name=self.name,output=fake_image,\
                    objects=starlist_name)

            starlist_ra,starlist_dec = self.wcs.wcs_pix2world(starlist_x,starlist_y,1)

            starlist_coords = coords.SkyCoord(ra=starlist_ra*u.degree,\
                dec=starlist_dec*u.degree)

            self.starlist_cat = se.seCatalog()
            self.starlist_cat.coords = starlist_coords
            self.starlist_cat.mags = starlist_mags
            self.starlist_cat.fwhm = starlist_fwhm
            self.starlist_cat.ellipticity=starlist_ellipticity

            self.runSE(fake_image,self.config_file,cat_name=recovered_cat_name)
            self.recovered_cat = se.seCatalog()
            self.recovered_cat.create_new_catalog_file(recovered_cat_name,filter='recovered')

            self.starlist_cat.mergeCatalog(self.recovered_cat,match_rad=match_rad)
            self.starlist_cat.mergeCatalog(self.phot,match_rad=match_rad)

            #create new photometry lists containing only the new stars
            old_phot = self.starlist_cat.mags['input'] > 99
            input_mags = self.starlist_cat.mags['input'][~old_phot]
            recovered_mags = self.starlist_cat.mags['recovered'][~old_phot]
            input_coords = self.starlist_cat.coords[~old_phot]
            input_pix_x,input_pix_y = self.wcs.wcs_world2pix(input_coords.ra*u.degree,input_coords.dec*u.degree,1)

            if overwrite or ini:
                key = str(i)
            else:
                key = str(max(map(int,self.fake_cats.keys()))+1)
            self.fake_cats[key] = np.vstack((input_mags,recovered_mags))
            self.fake_cats_x[key] = input_pix_x
            self.fake_cats_y[key] = input_pix_y

            un_rec = recovered_mags > 99
            self.unrec_cat= se.seCatalog()
            self.unrec_cat.coords = input_coords[un_rec]
            self.unrec_cat.id = np.zeros(un_rec.size)
            self.unrec_cat.mags = {'input':input_mags[un_rec]}

            self.rec_cat= se.seCatalog()
            self.rec_cat.coords = input_coords[~un_rec]
            self.rec_cat.id = np.zeros(un_rec.size)
            self.rec_cat.mags = {'input':input_mags[~un_rec]}


    def makePSF(self):
        self.runner.makePSF()

    def runSE(self,image_file,config_file,cat_name='sex.cat'):
        '''
        Run SE to get photometry
        '''
        sp.call('sex ' +image_file+' -CATALOG_NAME '+cat_name+' -c '+config_file,shell=True)


