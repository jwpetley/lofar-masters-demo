import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u  # type: ignore
from tqdm import tqdm
from glob import glob
import os

import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.visualization import (
    SqrtStretch,
    ImageNormalize,
    PercentileInterval,
)

from typing import Tuple


def open_table(filename, format="fits") -> Table:
    return Table.read(filename, format=format)


def clear_cutouts(path) -> None:
    fits_files = glob(path + "*.fits")
    png_files = glob(path + "*.png")

    for f in fits_files:
        os.remove(f)

    for f in png_files:
        os.remove(f)


def sort_table(t, column="Total_flux") -> Table:
    return t.sort(column, reverse=True)


def open_image(filename) -> Tuple[np.ndarray, WCS]:
    with fits.open(filename) as f:
        data: np.ndarray = f[0].data.squeeze().squeeze()  # Remove extra dimensions
        wcs: WCS = WCS(f[0].header).celestial  # Remove extra dimensions
        return data, wcs


def get_skycoord(source, ra_col="RA", dec_col="DEC") -> SkyCoord:
    ra = source[ra_col]
    dec = source[dec_col]
    return SkyCoord(ra=ra, dec=dec, unit="deg")


def process_table(t, img, wcs, n_sources=10) -> None:
    for source in tqdm(t[:n_sources]):
        c = get_skycoord(source)
        cutout = create_cutout(img, wcs, c)
        write_cutout(source, cutout)
        image_cutout(source, cutout)


def create_cutout(img, wcs, coord, size=5 * u.arcmin) -> Cutout2D:
    return Cutout2D(img, coord, size, wcs)


def write_cutout(source, cutout, dir="./CUTOUTS/") -> None:
    hdu = fits.PrimaryHDU(cutout.data, cutout.wcs.to_header())
    hdul = fits.HDUList([hdu])
    out_name = f"LoTSS_{source['RA']}_{source['DEC']}.fits"
    out_path = dir + out_name
    hdul.writeto(out_path, overwrite=True)


def image_cutout(source, cutout, dir="./CUTOUTS/") -> None:
    fig, ax = plt.subplots(subplot_kw=dict(projection=cutout.wcs))
    # Prepare the data
    norm = ImageNormalize(
        cutout.data, interval=PercentileInterval(99.9), stretch=SqrtStretch()
    )

    ax.imshow(cutout.data, origin="lower", norm=norm)

    out_name = f"LoTSS_{source['RA']}_{source['DEC']}.png"
    out_path = dir + out_name
    plt.savefig(out_path)
    plt.close()


def image_all_cutouts(t, dir="./CUTOUTS/") -> None:
    cutout_files = glob(dir + "*.fits")

    ncol = int(len(cutout_files) / 5)
    nrow = 5

    fig = plt.figure(figsize=(ncol + 1, nrow + 1))

    gs = gridspec.GridSpec(
        nrow,
        ncol,
        wspace=0.0,
        hspace=0.0,
        top=1.0 - 0.5 / (nrow + 1),
        bottom=0.5 / (nrow + 1),
        left=0.5 / (ncol + 1),
        right=1 - 0.5 / (ncol + 1),
    )

    transform = PercentileInterval(99.9) + SqrtStretch()

    for i, f in enumerate(cutout_files):
        ax = plt.subplot(gs[i])
        img, wcs = open_image(f)

        trans_data = transform(img)

        ax = plt.subplot(gs[i])
        ax.imshow(trans_data)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axis("off")

    plt.savefig("PLOTS/all_cutouts.png")


if __name__ == "__main__":
    table_path: str = "CATALOGUES/LOFAR_EDFN_v1.0.fits"
    image_path: str = (
        "IMAGES/image_full_ampphase_di_m.NS.int.restoredcorrectedflux.fits"
    )

    np.seterr(all="ignore")  # Disable warnings Annoying me

    clear_cutouts("./CUTOUTS/")

    edfn_lofar: Table = open_table(table_path)

    sort_table(edfn_lofar, column="Total_flux")

    img_data, wcs = open_image(image_path)

    process_table(edfn_lofar, img_data, wcs, n_sources=20)

    image_all_cutouts(edfn_lofar)
