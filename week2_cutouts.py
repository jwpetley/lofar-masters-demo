from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import (
    PercentileInterval,
    SqrtStretch,
    ZScaleInterval,
    LogStretch,
)
import os
from tqdm import tqdm


def findrms(mIn, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m = mIn[np.abs(mIn) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.0
    med = np.median(m)
    for i in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs((rms - rmsold) / rmsold) < diff:
            break
        rmsold = rms
    return rms


def load_image(file):
    with fits.open(file) as f:
        data = f[0].data.squeeze().squeeze()
        wcs = WCS(f[0].header).celestial
        return data, wcs


def download_high_res_EDFN(facet):
    """There are probably better ways to do this but oh well"""
    facet_fill = str(facet).zfill(2)
    url_05 = f"https://lofar-surveys.org/downloads/EDFN-HR/0.5res/images/poly_{facet_fill}_0.5asec_flux+astrometry_corr.fits"
    url_07 = f"https://lofar-surveys.org/downloads/EDFN-HR/0.7res/images/poly_{facet_fill}_int_corr.fits"

    out_05 = f"./IMAGES/EDFN_0.5_{facet_fill}.fits"
    out_07 = f"./IMAGES/EDFN_0.7_{facet_fill}.fits"

    user = "surveys"
    password = "150megahertz"

    os.system(f"wget -c --user {user} --password {password} {url_05} -O {out_05}")
    os.system(f"wget -c --user {user} --password {password} {url_07} -O {out_07}")


def cut_table_to_facet(table, facet_wcs):
    coords = SkyCoord(table["RA"], table["DEC"], unit="deg")
    inside_facet_mask = facet_wcs.footprint_contains(coords)

    sub_table = table[inside_facet_mask]
    print(f"{len(sub_table)} sources in facet")
    return sub_table


def generate_cutout(data, wcs, coord, size=3 * u.arcmin):
    cutout = Cutout2D(data, coord, size, wcs)
    return cutout


def plot_cutouts(cutouts, source):
    fig, axs = plt.subplots(2, 2)
    axs = axs.ravel()

    labels = ["0.5 arcsec", "0.7 arcsec", "1.5 arcsec", "6 arcsec"]

    transform = LogStretch()

    for i, cutout in enumerate(cutouts):
        rms = findrms(cutout.data)

        trans_data = transform(cutout.data)

        axs[i].imshow(trans_data, origin="lower", vmin=4 * rms, vmax=100 * rms)
        axs[i].axis("off")
        axs[i].set_title(labels[i])

    plt.suptitle(
        f"Source id: {source['Source_id']}\nTotal flux: {source['Total_flux']:.2f} Jy"
    )
    plt.tight_layout()
    plt.show()


def check_valid_cutout(cutout):
    nan_count = np.count_nonzero(np.isnan(cutout.data))

    if nan_count > 0.8 * np.size(cutout.data):
        return False
    else:
        return True


if __name__ == "__main__":
    facet = 0
    facet_fill = str(facet).zfill(2)

    file_05 = f"./IMAGES/EDFN_0.5_{facet_fill}.fits"
    file_07 = f"./IMAGES/EDFN_0.7_{facet_fill}.fits"
    file_15 = "./IMAGES/4days_1.5asec_flux+astro-corr.fits"
    file_60 = "./IMAGES/image_full_ampphase_di_m.NS.int.restoredcorrectedflux.fits"

    if not os.path.isfile(file_05):
        download_high_res_EDFN(facet)

    if not os.path.isfile(file_07):
        download_high_res_EDFN(facet)

    edfn_table = Table.read("CATALOGUES/LOFAR_EDFN_v1.0.fits")

    data_05, wcs_05 = load_image(file_05)
    data_07, wcs_07 = load_image(file_07)
    data_15, wcs_15 = load_image(file_15)
    data_60, wcs_60 = load_image(file_60)

    edfn_table_facet = cut_table_to_facet(edfn_table, wcs_07)

    edfn_table_facet.sort("Total_flux", reverse=True)

    for source in tqdm(edfn_table_facet[0:10]):
        coord = SkyCoord(source["RA"], source["DEC"], unit="deg")

        cutout_05 = generate_cutout(data_05, wcs_05, coord)
        cutout_07 = generate_cutout(data_07, wcs_07, coord)
        cutout_15 = generate_cutout(data_15, wcs_15, coord)
        cutout_60 = generate_cutout(data_60, wcs_60, coord)

        if check_valid_cutout(cutout_05):
            plot_cutouts([cutout_05, cutout_07, cutout_15, cutout_60], source)
