#first import modules
from astropy.io import fits
from numpy import array, float32, float16, radians, resize, linspace, digitize
from numpy import sqrt, degrees, arctan, mean, median, std, argsort, transpose, savetxt, arange, tile, where,  savez, load
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import scipy.signal
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import SmoothBivariateSpline
import warnings


def mapwavelength(trace, wavemap, mode='spline2d'):
    """
    Compute the wavelength along the center of the trace, to be run after
    the HeNeAr_fit routine.

    Parameters
    ----------
    trace : 1-d array
        The spatial positions (Y axis) corresponding to the center of the
        trace for every wavelength (X axis), as returned from ap_trace
    wavemap : bivariate spline object or image-like wavelength map
        The wavelength evaluated at every pixel, output from HeNeAr_fit
        Type depends on mode parameter.
    mode : str, optional
        Which mode was used to generate the 2D wavelength solution in
        HeNeAr_fit(), and specifically in lines_to_surface()?
        Options include: poly, spline, spline2d (Default is 'spline2d')

    Returns
    -------
    trace_wave : 1d array
        The wavelength vector evaluated at each position along the trace
    """
    # use the wavemap from the HeNeAr_fit routine to determine the wavelength along the trace
    if mode=='spline2d':
        trace_wave = wavemap.ev(np.arange(len(trace)), trace)

    elif mode=='poly' or mode=='spline':
        trace_wave = np.zeros_like(trace)
        for i in range(len(trace)):
            trace_wave[i] = np.interp(trace[i], range(wavemap.shape[0]), wavemap[:,i])
    ## using 2d polyfit
    # trace_wave = polyval2d(np.arange(len(trace)), trace, wavemap)
    return trace_wave

def AirmassCor(obj_wave, obj_flux, airmass, airmass_file='atmoexan.dat'):
    """
    Correct the spectrum based on the airmass

    Parameters
    ----------
    obj_wave : 1-d array
        The 1-d wavelength array of the spectrum
    obj_flux : 1-d or 2-d array
        The 1-d or 2-d flux array of the spectrum
    airmass : float
        The value of the airmass, not the header keyword.
    airmass_file : str, {'apoextinct.dat', 'ctioextinct.dat', 'kpnoextinct.dat', 'ormextinct.dat'}
        The name of the airmass extinction file. This routine assumes
        the file is stored in the resources/extinction/ subdirectory.
        Available files are (Default is apoextinct.dat)

    Returns
    -------
    The flux array
    """
    # read in the airmass extinction curve
   # extinction_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
   #                           'resources/extinction')
   # if len(airmass_file)==0:
   #     air_wave, air_cor = np.genfromtxt(os.path.join(extinction_dir, airmass_file),
   #                                    unpack=True,skip_header=2)
   # else:
        #print('> Loading airmass library file: '+airmass_file)
        # print('  Note: first 2 rows are skipped, assuming header')
    air_wave, air_cor = np.genfromtxt(airmass_file, unpack=True,skip_header=1)
    # air_cor in units of mag/airmass

    airmass_ext = 10.0**(0.4 * airmass *np.interp(obj_wave, air_wave, air_cor))
    # arimas_ext is broadcast to obj_flux if it is a 2-d array
    return obj_flux * airmass_ext

def _WriteSpec(spec, wfinal, file_name):
    # write the final spectrum out
    fout = open(file_name,'w')
    fout.write('#  This file contains the final extracted (wavelength, spec normalised) data \n')
    for k in range(len(wfinal)):
        fout.write(str(wfinal[k]) + '  ' + str(spec[k]) + '\n')
    fout.close()
    #print(file_name+'_atmocorr.txt saved !')
    return


# fichier fits contenant la table de calibration en longueur d'onde (obtenue via arcline_calib.py)
hdu = fits.open('wavefit2_glo.fits')
wfit=hdu[0].data
hdu.close(closed=True)
print(wfit.shape)


#On ouvre le fichier .fits correspondant au spectre à corriger
filename = "AF708430p0.9_mono"
spectra_fits = fits.open(filename + ".fits")

#On trace ce spectre et la masse d'air
spectra_hdr = spectra_fits[0].header
spectra = spectra_fits[0].data
airmass = spectra_fits[0].header['AIRMASS']
print("This is the spectrum airmass")
print(airmass)
wfinal = mapwavelength(np.arange(spectra.shape[0]), wfit[11:-15,:], mode='poly')   # calibration en lbd
print(wfinal.shape)





#On utilise la fonction AirmassCor pour corriger l'atmosphère
#On trace les spectres avant / après correction atmosphérique


    # Génère l'axe des longueurs d'onde la correction à apporter
spectre_corrige=AirmassCor(wfinal,spectra,airmass,airmass_file=('ekar_sloanfilter_ext'+'.txt')) #qui associe à chaque longueur d'onde
plt.figure()
plt.plot(wfinal,spectra,label="spectre non corrigé")
plt.legend()    
plt.xlabel("longueur d'onde")
plt.ylabel('flux non corrigé')
plt.show()

plt.figure()
plt.plot(wfinal,spectre_corrige,label="spectre corrigé")
plt.legend()    
plt.xlabel("longueur d'onde")
plt.ylabel('flux corrigé')
plt.show()

_WriteSpec(spectre_corrige, wfinal, filename + '_aircorr.txt') #remplacez ext_spec_filt par le nom du spectre corrigé

ext_spec_filt = spectre_corrige
        # -- Normalisation par la moyenne du spectre entre 6380 et 6420 A 
arg_wfinal_norm_left = np.searchsorted(wfinal, 6380, side='left')    # recherche des indices correspondant a lbd = 5480 A
arg_wfinal_norm_right = np.searchsorted(wfinal,6420,side='right')
ext_spec_filt_norm = ext_spec_filt / np.mean(ext_spec_filt[arg_wfinal_norm_left:arg_wfinal_norm_right])   # normalisation du spectre par la moyenne de l'amplitude du spectre sur 5480 < lbd < 5520 A
        # spec_filt_norm = spec_filt / np.mean(spec_filt[arg_wfinal_norm:arg_wfinal_norm+21])

# On trace le spectre corrigé de la masse d'air et normalisé
plt.figure()
plt.plot(wfinal, ext_spec_filt_norm, label = 'with airmass correction + normalisation')
plt.legend()
plt.xlabel('Wavelength (angstrom)')
plt.ylabel('Intensity (ADU)')
plt.title("Spectre calibré et corrigé de la masse d'air atmosphérique")
plt.show()

# On enregistre le spectre calibré spectralement, corrigé de la masse d'air et normalisé
_WriteSpec(ext_spec_filt_norm, wfinal, filename + '_aircorr_norm.txt')
