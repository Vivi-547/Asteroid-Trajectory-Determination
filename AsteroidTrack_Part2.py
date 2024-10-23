from Master_OD_Functions import *
import numpy as np
from odlib.odlib import *
from astropy.io import fits
import matplotlib.pyplot as plt

# Given Quantities
k = 0.0172020989484 #Gaussian gravitational constant
cAU = 173.144643267 #speed of light in au/(mean solar)day
eps = math.radians(23.4374) #Earth's obliquity

# MONTECARLO - add RMS-tweaked RA and DEC to each observation array
ras1 = []
ras2 = []
ras3 = []

decs1 = []
decs2 = []
decs3 = []
for i in range(3):
    if (i == 0):
        table = fits.open('corr1.fits')[1].data
    elif (i == 1):
        table = fits.open('corr2.fits')[1].data
    else:
        table = fits.open('corr3.fits')[1].data

    field_ra = table.field_ra
    field_dec = table.field_dec
    index_ra = table.index_ra
    index_dec = table.index_dec

    mean_ra = []
    mean_dec = []
    for i in range(int(table.shape[0])):
        mean_ra.append((field_ra[i] - index_ra[i])**2)
        mean_dec.append((field_dec[i] - index_dec[i])**2)

    rms_ra = math.sqrt(np.mean(mean_ra))
    rms_dec = math.sqrt(np.mean(mean_dec))

posit, veloc, semi_major_arr, eccent_arr, inclin_arr, longit_asc_arr, arg_perih_arr, mean_anom_arr = Final("viviInput.txt", rms_ra, rms_dec)
print("real final posit, veloc", posit, veloc)

# Create matrix [array, mean, std] for each orbital element
# Order is semi_major(A), eccent(E), inclin(I), longit asc(O), arg perih(o), mean anom(M)
orbitalelem = [
    [],
    [],
    [],
    [],
    [],
    []
]
orbitalelem[0] = [semi_major_arr, np.mean(semi_major_arr), np.std(semi_major_arr)]
orbitalelem[1] = [eccent_arr, np.mean(eccent_arr), np.std(eccent_arr)]
orbitalelem[2] = [inclin_arr, np.mean(inclin_arr), np.std(inclin_arr)]
orbitalelem[3] = [longit_asc_arr, np.mean(longit_asc_arr), np.std(longit_asc_arr)]
orbitalelem[4] = [arg_perih_arr, np.mean(arg_perih_arr), np.std(arg_perih_arr)]
orbitalelem[5] = [mean_anom_arr, np.mean(mean_anom_arr), np.std(mean_anom_arr)]

predictedVal = [3.423119774454217E+00, 6.770088849489315E-01, 5.616390117735835E+01, 2.530523773999922E+02, 3.008338981468623E+02, 1.303659573706276E+01]

# print(f"semi major error = {error(a, 3.423119774454217E+00)}")
#         print(f"eccentricity error = {error(e, 6.770088849489315E-01)}")
#         print(f"inclination error = {error(incline, 5.616390117735835E+01)}")
#         print(f"logitude of ascending node error = {error(Big_Omega, 2.530523773999922E+02)}")
#         print(f"True anomaly error = {error(Nu, 7.243923573300094E+01)}")
#         print(f"argument of perihelion error = {error(small_Omega, 3.008338981468623E+02)}")
#         print(f"Mean Anomaly error = {error(M,1.303659573706276E+01)}")

# Histogram
#fig, ax = plt.subplots(6,1)
label = ['Semi-Major Axis (a)', 'Eccentricity (e)', 'Inclination (i)', 'Longitude of Ascending Node ('r"$\Omega$"')', 'Argument of Perihelion ('r"$\omega$"')', 'Mean Anomaly (M)']
#range_arr = [np.linspace(0,5,50), np.linspace(0,1,50), np.linspace(50,60,50), np.linspace(240,260,50), np.linspace(290,310,50), np.linspace(10,15,50)]
range_arr = [np.linspace(2,5,50), np.linspace(0,3,50), np.linspace(55,58,50), np.linspace(251,254,50), np.linspace(299,302,50), np.linspace(12,15,50)]
for i in range(6):
    fig, ax = plt.subplots(1,1)
    mean = orbitalelem[i][1]
    print("Mean " + str(label[i]) + ": " + str(mean))
    std = orbitalelem[i][2]
    print("Standard Deviation " + str(label[i]) + ": " + str(mean))
    omega1 = orbitalelem[i][0]
    def gaussian(x, a, mean, std):
        return 1./(np.sqrt(np.pi * 2) * std)*np.exp(-(x-mean)**2/(2*std**2))
    
    ax.hist(omega1, alpha=0.75, bins = 25, zorder=2, histtype="step", linewidth=3, density=True)
    ax.axvline(predictedVal[i], color="black", linestyle="--", label="JPL Value= " + str(predictedVal[i]))
    span = np.linspace(mean - std, mean + std, 50)  # the range of mean +/- std
    xrange = range_arr[i] # the full xspan of the plot

    ax.plot(xrange, gaussian(xrange, None, mean, std), color="black")
    ax.fill_between(span, np.zeros(50), gaussian(span, None, mean, std),
                    alpha=0.50, color="orange", zorder=3, label = r"$\sigma$= %.7F" % std)

    ax.set_xlabel(label[i], fontname="Times New Roman", fontsize=16)

    # Legend just displays the labels
    ax.legend()

plt.show()
