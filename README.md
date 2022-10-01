# fconvo
`fconvo` is a Python module to obtain accurate synthetic photometry given input spectra and filter curves, taking into account whether the filters are photon-counting or energy-counting. It was written by Nadia Zakamska. This module is written for Python 3, and it depends only on `numpy`, though the tutorial also uses `astropy, matplotlib, scipy`. If you use `fconvo` in your research, please cite https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4191Z/abstract.

There are three functions in `fconvo`: 
1. `rebin` rebins the spectrum onto another wavelength grid using linear interpolation. Its functionality is similar to many others, including `interp1d` from `scipy`. Unlike `interp1d`, `rebin` does not break down if the wavelength coverage of the original spectrum is inadequate. Although a flag coverage is returned, all warnings are suppressed, and it's the user's responsibility to check the coverage of the spectral templates used in synthetic photometry. 
2. `fnu_ave` and `flam_ave` average the flux taking into account filter transmission functions and using equations appropriate for photon-counting and energy-counting instruments. 

Among other examples, `fconvo` uses the SWIRE library of galaxy templates (http://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html) which are included with `fconvo` for ease of use. If you download these templates here and use them in your research, please cite the original SWIRE papers (e.g., https://ui.adsabs.harvard.edu/abs/2007ApJ...663...81P/abstract and references therein). 

Equations used for photon-counting and energy-counting filters are from here: https://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/. In case this webpage disappears, we also provide a pdf printout of this webpage as `intro_filter_averaging.pdf`. 

A limited set of filter curves is supplied with this distribution in the FIL_CURVE directory. The filter curves have their provenance and type (photon-counting vs energy-counting) listed in the headers. 

We provide a tutorial Jupyter notebook for this module, `tutorial.ipynb`. The tutorial consists of preliminaries (imports) and four examples. 

In Example 1, a spectrum is rebinned onto a new wavelength grid and shown using `matplotlib` with the original spectrum and two other implementations of rebinning from SciPy. 

In Example 2, AB colors of a galaxy placed at different redshifts are computed. 

In Example 3, average F_lambda is shown together with the original spectrum. 

In Example 4, synthetic AB photometry of a white dwarf is computed and compared with observed SDSS photometry. In this example, the spectral coverage of one of the bands is inadequate, resulting in a disagreement between synthetic and observed photometry. 

## Contact

Feel free to contact me: zakamska@jhu.edu

## References

The module was originally developed for stitching the orders of Spitzer spectra to match WISE photometry (Hill and Zakamska, 2014, https://ui.adsabs.harvard.edu/abs/2014MNRAS.439.2701H/abstract) and then used for synthetic photometry at other wavelengths (Zakamska et al. 2016, https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4191Z/abstract)

