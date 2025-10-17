# LOFAR Master's Code Demos


## Week 1

- `week1_spectral_index.py` requires LoTSS, VLASS and NVSS catalogues
- `week1_cutouts.py` requires low-res EDFN image and catalogue - [https://lofar-surveys.org/deepfields_public_edfn.html](https://lofar-surveys.org/deepfields_public_edfn.html)

To make the file structure work you need to have CATALOGUES, CUTOUTS, IMAGES and PLOTS directories.


You also need to have `astropy`, `numpy`, `scipy` and `matplotlib` installed

e.g `pip install astropy numpy scipy matplotlib`

You can run the code by  :

```bash
python week1_spectral_index.py # Cross match catalogues and plot spectral indices

python week1_cutouts.py # Make cutouts of the brightest sources in EDFN
```


The second script will produce an image like this:

![LOFAR Cutouts](./all_cutouts.png)