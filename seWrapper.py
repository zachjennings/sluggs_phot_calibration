import numpy as np
import subprocess as sp
import os
import pickle
import astropy.coordinates as coords
import astropy.units as u
from sklearn.svm import SVC as svc
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.grid_search import GridSearchCV
from astropy.wcs import WCS
from astropy.io import fits
from pyraf.iraf import artdata

class seCatalog(object):
    def __init__(self):
        self.mags = dict()
        self.ellipticity = dict()
        self.fwhm = dict()
        self.flags = dict()
        self.merr = dict()
        
    def readSE(self,output_file):
        '''
        Read in the given output file
        '''
        new_catalog = np.loadtext(output_file)

    def create_new_catalog_file(self,output_file,filter='g'):
        new_catalog = np.loadtxt(output_file,comments='#')
        bad = new_catalog[:,2] > 90
        new_catalog=new_catalog[~bad,:]
        ra = new_catalog[:,0]
        dec = new_catalog[:,1]
        c = coords.SkyCoord(ra = ra*u.degree,dec = dec*u.degree)
        self.coords = c
        self.mags[filter] = new_catalog[:,2]
        self.fwhm[filter] = new_catalog[:,3]
        self.ellipticity[filter] = new_catalog[:,4]
        self.flags[filter] = new_catalog[:,5]
        
    def create_new_catalog_arrs(self,ra,dec,fil,fwhm=None,mags=None,merr=None,flags=None,ellipticity=None):
        """
        Create a new catalog. Ra and Dec are n x 1 np arrays. 
        """
        if mags is None:
            mags = np.zeros(ra.size)
        if fwhm is None:
            fwhm = np.zeros(ra.size)
        if merr is None:
            merr = np.zeros(ra.size)
        if flags is None:
            flags = np.zeros(ra.size)
        if ellipticity is None:
            ellip = np.zeros(ra.size)
            
        c = coords.SkyCoord(ra = ra*u.degree,dec = dec*u.degree)
        self.coords = c
        self.mags[fil] = mags
        self.fwhm[fil] = fwhm
        self.merr[fil] = merr
        self.ellipticity[fil] = ellipticity
        self.flags[fil] = flags

    def mergeCatalog(self,new_catalog,match_rad=0.5*u.arcsec):
        '''
        Merge a given seCatalog into the current one.
        '''

        full_old,full_new,d2d,d3d = new_catalog.coords.search_around_sky(self.coords,match_rad)
        #logic to find the minimum seperation for a given "new"
        old = np.array([],dtype=int)
        new = np.array([],dtype=int)

        #Sort all arrays by d2d, so that we know the first value will be the correct one
        d_sort = np.argsort(d2d)
        full_old = full_old[d_sort]
        full_new = full_new[d_sort]
        d2d = d2d[d_sort]

        self.full_old = full_old
        self.full_new = full_new
        self.d2d = d2d

        for i in np.arange(d2d.size):
            this_new = full_new[i]
            this_old = full_old[i]

            #if either match has already been "claimed", reject it
            if this_new not in new and this_old not in old:
                #old = np.append(old,this_old[np.argmin(d2d.arcsec[full_new == i])])
                #new = np.append(new,closest_new)
                old = np.append(old,this_old)
                new = np.append(new,this_new)

#       for i in new_catalog.coords:
#           min_dist = np.min(i.seperation(self.coords))
#           if min_dist < match_rad:
#               this_new = np.array([min_dist])
#
#           np.append(old,this_old)
#           np.append(new,this_new[np.argmin(d2d)])

        #logic to find nearest matches for given coordinates
#       sort = np.argsort(new)
        self.new=new
        self.old = old
        arr_old = np.arange(self.coords.ra.size,dtype=int)
        self.arr_old = arr_old
        arr_new = np.arange(new_catalog.coords.ra.size,dtype=int)
        self.arr_new = arr_new
        idx_old = np.array([i in old for i in arr_old])
        idx_new = np.array([i in new for i in arr_new])
        if np.sum(idx_old) < 1:
            old = np.zeros(self.coords.ra.size,dtype=bool)
            idx_old = np.zeros(self.coords.ra.size,dtype=bool)
            new = np.zeros(new_catalog.coords.ra.size,dtype=bool)
            idx_new = np.zeros(new_catalog.coords.ra.size,dtype=bool)
            print('No matches found')
        #idx_old = np.in1d(np.arange(self.coords.ra.size,dtype=int),old)
        #idx_new = np.in1d(np.arange(new_catalog.coords.ra.size,dtype=int),new)
        self.idx_old = idx_old
        self.idx_new = idx_new
        #select out the matched coordinates

        unmatched_old_coords = self.coords[~idx_old]
        unmatched_new_coords = new_catalog.coords[~idx_new]

        if np.any(idx_old):
            matched_coords = self.coords[old]
            new_ra = np.concatenate([matched_coords.ra.deg,unmatched_old_coords.ra.deg,unmatched_new_coords.ra.deg])
            new_dec = np.concatenate([matched_coords.dec.deg,unmatched_old_coords.dec.deg,unmatched_new_coords.dec.deg])
            self.matched_coords = matched_coords

        else:
            new_ra = np.concatenate([unmatched_old_coords.ra.deg,unmatched_new_coords.ra.deg])
            new_dec = np.concatenate([unmatched_old_coords.dec.deg,unmatched_new_coords.dec.deg])

        self.unmatched_old_coords = unmatched_old_coords
        self.unmatched_new_coords = unmatched_new_coords
        #m_ind = matched_coords.ra.size
        #um_ind_old = m_len + unmatched_old_coords.ra.size
        #um_ind_new = m_len + um_ind_old + unmatched_new_coords.ra.size
        self.new_ra = new_ra
        self.new_dec = new_dec
        #first concetenate the RAs and DECs of the two catalogs
        self.old_coords = self.coords
        self.coords = coords.SkyCoord(ra = new_ra*u.degree,dec=new_dec*u.degree)
        if np.any(idx_old):
            detected = np.ones(matched_coords.ra.size,dtype=bool)
        else:
            detected=np.array([])
        undetected = np.zeros(unmatched_old_coords.ra.size + unmatched_new_coords.ra.size,dtype=bool)
        self.all_detected = np.concatenate([detected,undetected])

        #Magnitudes:
        #iterate over the keys for each catalog in magnitude
        for i in self.mags.keys():
            matched_mags = self.mags[i][old]
            unmatched_mags = self.mags[i][~idx_old]
            non_detections = np.zeros(unmatched_new_coords.ra.size) + 99.999
            self.non_detections_new = non_detections
            self.mags[i] = np.concatenate([matched_mags,unmatched_mags,non_detections])

        for i in new_catalog.mags.keys():
            matched_mags = new_catalog.mags[i][new]
            unmatched_mags = new_catalog.mags[i][~idx_new]
            non_detections = np.zeros(unmatched_old_coords.ra.size) + 99.999
            self.non_detections_old = non_detections
            self.mags[i] = np.concatenate([matched_mags,non_detections,unmatched_mags])

        #FWHM:
        #next iterate over the keys for each catalog in FWHM
        for i in self.fwhm.keys():
            matched_fwhm = self.fwhm[i][old]
            unmatched_fwhm = self.fwhm[i][~idx_old]
            non_detections = np.zeros(unmatched_new_coords.ra.size) + 99.999
            self.fwhm[i] = np.concatenate([matched_fwhm,unmatched_fwhm,non_detections])

        for i in new_catalog.mags.keys():
            matched_fwhm = new_catalog.fwhm[i][new]
            unmatched_fwhm = new_catalog.fwhm[i][~idx_new]
            non_detections = np.zeros(unmatched_old_coords.ra.size) + 99.999
            self.fwhm[i] = np.concatenate([matched_fwhm,non_detections,unmatched_fwhm])

        #Ellipticity:
        #next iterate over the keys for each catalog in ellipticity
        for i in self.ellipticity.keys():
            matched_ellip = self.ellipticity[i][old]
            unmatched_ellip = self.ellipticity[i][~idx_old]
            non_detections = np.zeros(unmatched_new_coords.ra.size) + 99.999
            self.ellipticity[i] = np.concatenate([matched_ellip,unmatched_ellip,non_detections])

        for i in new_catalog.mags.keys():
            matched_ellip = new_catalog.ellipticity[i][new]
            unmatched_ellip = new_catalog.ellipticity[i][~idx_new]
            non_detections = np.zeros(unmatched_old_coords.ra.size) + 99.999
            self.ellipticity[i] = np.concatenate([matched_ellip,non_detections,unmatched_ellip])
            
        for i in self.flags.keys():
            matched_flags = self.flags[i][old]
            unmatched_flags = self.flags[i][~idx_old]
            non_detections = np.zeros(unmatched_new_coords.ra.size) + 99.999
            self.flags[i] = np.concatenate([matched_flags,unmatched_flags,non_detections])
            
        for i in new_catalog.flags.keys():
            matched_flags = new_catalog.flags[i][new]
            unmatched_flags = new_catalog.flags[i][~idx_new]
            non_detections = np.zeros(unmatched_old_coords.ra.size) + 99.999
            self.flags[i] = np.concatenate([matched_flags,non_detections,unmatched_flags])

        truth = np.zeros((self.coords.ra.size,len(self.mags.keys())),dtype=bool)
        keys = self.mags.keys()
        for i in range(len(list(keys))):
            truth[:,i] = self.mags[list(keys)[i]] < 99

        self.all_detected = np.all(truth,axis=1)
        self.id = np.arange(self.all_detected.size)

    def getXY(self,image):
        '''
        Return coordinates of sources in pixel coordinates, if desired.

        Inputs:
            Fits: fits image from which to take WCS

        Returns:
            x,y : x and y pixel coordinates for each source
        '''
        wcs = WCS(image)
        return wcs.wcs_world2pix(self.coords.ra*u.degree,self.coords.dec*u.degree,1)

    def coordPrint(self,fileroot,size='1.0"',color='green',write_comment=True,select='',det=True):
        '''
        select = only pick detected mags filter
        det = selected detected or non-detected
        '''
        with open(fileroot+'_matched_only.reg',"w+") as the_file:
            the_file.write('fk5 \n')
            if select and det:
                good = np.where(self.mags[select] < 99)[0]
            elif not det:
                good = np.where(self.mags[select] > 99)[0]
            else:
                good = np.ones(self.coords.ra.size,dtype='bool')
            for i in np.arange(good.size):
                #determine the text of the comment if specified
                if write_comment:
                    keys = self.mags.keys()
                    if select:
                        the_key = select
                    else:
                        the_key = keys[1]

                    mag = self.mags[the_key][good[i]]

                    comment='# text = {ID='+str(self.id[good[i]])+' '+\
                    str(the_key)+'='+str(mag)+'}'
                else:
                    comment = ''
                the_file.write('circle '+\
                str(self.coords.ra.deg[good[i]])+' '+str(self.coords.dec.deg[good[i]])+' '+ \
                size + comment + '\n')



class fakeStarTests(object):
    def __init__(self,image,phot_catalog,pix_scale=0.2,gain=1.0,rd_noise=0.0,\
        zpt=30.,seeing=3.5,config_file='foo.sex',write_true_phot=False):
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
                        magzero=self.zpt,gain=self.gain,rdnoise=self.rd_noise,radius=self.seeing,background=0.0)


    def makeStarList(self,name,n_stars=1000,min_mag=24,max_mag=28,x_range=(0,1000),y_range=(0,1000),write_true_phot=False):
        '''
        Use pyraf artdata.starlist to make an input starlist
        '''
        #artdata.starlist(name,nstars=n_stars,xmax=self.header['naxis1'],ymax=self.header['naxis2'],\
        #    minmag=min_mag,maxmag=max_mag,luminosity='uniform',lseed='INDEF')

        mags = np.random.uniform(low=min_mag,high=max_mag,size=n_stars)
        x = np.random.uniform(high=self.header['naxis1'],size=n_stars)
        y = np.random.uniform(high=self.header['naxis2'],size=n_stars)

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

    def genFakeCatalogs(self,n_iter=1,min_mag=22.,max_mag=26.,n_stars=1000,overwrite=False):
        '''
        Perform photometry on the fake catalogs
        '''
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


            starlist_x,starlist_y,starlist_mags = self.makeStarList(starlist_name,min_mag=min_mag,max_mag=max_mag,n_stars=n_stars)
            starlist_mags={'input':starlist_mags}
            starlist_fwhm = {'input':np.zeros(starlist_x.size)+99.999}
            starlist_ellipticity = {'input':np.zeros(starlist_x.size)+99.999}

            self.placeStars(name=self.name,output=fake_image,\
                    objects=starlist_name)

            starlist_ra,starlist_dec = self.wcs.wcs_pix2world(starlist_x,starlist_y,1)

            starlist_coords = coords.SkyCoord(ra=starlist_ra*u.degree,\
                dec=starlist_dec*u.degree)

            self.starlist_cat = seCatalog()
            self.starlist_cat.coords = starlist_coords
            self.starlist_cat.mags = starlist_mags
            self.starlist_cat.fwhm = starlist_fwhm
            self.starlist_cat.ellipticity=starlist_ellipticity

            self.runSE(fake_image,self.config_file,cat_name=recovered_cat_name)
            self.recovered_cat = seCatalog()
            self.recovered_cat.create_new_catalog_file(recovered_cat_name,filter='recovered')

            self.starlist_cat.mergeCatalog(self.recovered_cat)
            self.starlist_cat.mergeCatalog(self.phot)

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
            self.unrec_cat=seCatalog()
            self.unrec_cat.coords = input_coords[un_rec]
            self.unrec_cat.id = np.zeros(un_rec.size)
            self.unrec_cat.mags = {'input':input_mags[un_rec]}

            self.rec_cat=seCatalog()
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






