# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 09:03:47 2025

@author: lenam
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import math
        
def _WriteSpec(spec, wfinal, file_name):
    """
    Writes the final spectrum out.
    """
    fout = open(file_name,'w')
    fout.write('#  This file contains the final extracted (wavelength, spec normalised, airmass, shifted, divised by star, filtered, cut) data \n')
    for k in range(len(wfinal)):
        fout.write(str(wfinal[k]) + '  ' + str(spec[k]) + '\n')
    fout.close()
    #print(file_name+'_atmocorr.txt saved !')
    return

#On ouvre les spectres de l'astéroïde et de l'analogue solaire, avec la correction atmosphérique appliquée
data_lex = np.loadtxt('AF708430p0.9_mono_aircorr_norm.txt')
data_star = np.loadtxt('AF708482_mono_aircorr_norm.txt')

#On sépare la longueur d'onde du flux
wavelenght_lex = (data_lex[:,0])
data_lex_y = (data_lex[:,1])
data_star_y = (data_star[:,1])

#On divise par l'analogue solaire
spectre_corrige = data_lex_y/data_star_y

#On trace le spectre corrigé
plt.figure()
plt.plot(wavelenght_lex,spectre_corrige,label="spectre_corrige")
plt.xlabel("Wavelenght (A)")
plt.ylabel("Intensity (ADU)")
plt.title("Spectre divisé par l'analogue solaire")
plt.show()

#On shift de quelques pixels
star_y = data_star_y
x = wavelenght_lex
y = data_lex_y

starysave = star_y
xsave = x
ysave = y

shift =3 #modification de la valeur du shift

star_y=starysave[:-shift] #on applique le shift à l'étoile
x=xsave[shift:]
y=ysave[shift:]

y_shifted = y/star_y #on divise à nouveau l'objet par l'analogue solaire

#On trace la comparaison avant/après correction du shift
plt.figure(figsize=(7,5))
plt.plot(x,y_shifted, label = '(2004)Lexell/SA102-1081 with shift correction', color='black')
plt.plot(wavelenght_lex, spectre_corrige, label = '(2004)Lexell/SA102-1081 without shift correction', color='grey', alpha=0.4)
plt.xlabel('Wavelegth (A)')
plt.ylabel('Intensity (ADU)')
plt.title("Comparaison avant/après correction du shift")
plt.legend()
plt.show

#On applique le filtre médian
window_size = 7 # Doit être un nombre impair
y_fin = medfilt(y_shifted, kernel_size=window_size)

#On coupe les premières et dernières longueurs d'ondes
cut = math.ceil((10251-9200)/5.4)
x_fin2 = x[:-cut]
y_fin2 = y_fin[:-cut]

#On convertie les ångström en micromètre
for k in range(len(x_fin2)):
    x_fin2[k] = x_fin2[k]*(10**(-4))
    
#On trace de le spectre divisé, filtré et coupé
plt.figure()
plt.plot(x_fin2,y_fin2,label="spectre filtré et coupé")
plt.xlabel("Wavelenght (µm)")
plt.ylabel("Intensity (ADU)")
plt.title("Spectre final filtré, coupé et divisé par l'analogue solaire")
plt.show()

_WriteSpec(y_fin2, x_fin2,'AF708430p0.9_mono_aircorr_norm_filter_cut.txt') #modifiez avec le nom que vous souhaitez donner au spectre