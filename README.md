# interspec
C code (using Pgplot interactively) for the analysis and visualisation of absorption spectra from radio telescopes

This will read in a radio-band spectrum as an asci file and recalibrate the axes to velocity or frequency (abscissa) and flux density/optical depth (ordinate). The spectrum can be resampled to any arbitrary velocity resolution and a polynomial fit subtracted from the bandpass. Once smoothed and bandpass subtracted, the rms noise level, peak optical depth, integrated optical depth, integrated flux and velocity (redshift) of the spectum will be given. These values are more physically justified than the usual practice of assuming a number of Gaussian fits (although these can also be performed).

There is also an option to produce a postscript plot.

![](https://raw.githubusercontent.com/steviecurran/interspec/refs/heads/main/1555-140.dat-2_poly_2.png)
