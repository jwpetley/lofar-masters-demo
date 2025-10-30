from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
from astropy.nddata import Cutout2D
from tqdm import tqdm
import os
from astropy.visualization import (
    PercentileInterval,
    SqrtStretch,
    ZScaleInterval,
    LogStretch,
)
from io import StringIO
import requests


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


# Function taken from https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812251/PS1+Image+Cutout+Service#PS1ImageCutoutService-PythonExampleScript
def getimages(tra, tdec, size=240, filters="grizy", format="fits", imagetypes="stack"):
    """Query ps1filenames.py service for multiple positions to get a list of images
    This adds a url column to the table to retrieve the cutout.

    tra, tdec = list of positions in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    format = data format (options are "fits", "jpg", or "png")
    imagetypes = list of any of the acceptable image types.  Default is stack;
        other common choices include warp (single-epoch images), stack.wt (weight image),
        stack.mask, stack.exp (exposure time), stack.num (number of exposures),
        warp.wt, and warp.mask.  This parameter can be a list of strings or a
        comma-separated string.

    Returns an astropy table with the results
    """

    ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"

    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    # if imagetypes is a list, convert to a comma-separated string
    if not isinstance(imagetypes, str):
        imagetypes = ",".join(imagetypes)
    # put the positions in an in-memory file object
    cbuf = StringIO()
    cbuf.write("\n".join(["{} {}".format(ra, dec) for (ra, dec) in zip(tra, tdec)]))
    cbuf.seek(0)
    # use requests.post to pass in positions as a file
    r = requests.post(
        ps1filename, data=dict(filters=filters, type=imagetypes), files=dict(file=cbuf)
    )
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")

    urlbase = "{}?size={}&format={}".format(fitscut, size, format)
    tab["url"] = [
        "{}&ra={}&dec={}&red={}".format(urlbase, ra, dec, filename)
        for (filename, ra, dec) in zip(tab["filename"], tab["ra"], tab["dec"])
    ]
    return tab


def download_cutouts(ps1_table, outdir="./CUTOUTS/"):
    for row in tqdm(ps1_table):
        ra = row["RA"]
        dec = row["DEC"]
        url = row["url"]

        outname = outdir + f"ps1_{ra}_{dec}.fits"

        if not os.path.isfile(outname):
            r = requests.get(url)
            open(outname, "wb").write(r.content)


def write_cutout(source, cutout, dir="./CUTOUTS/") -> None:
    hdu = fits.PrimaryHDU(cutout.data, cutout.wcs.to_header())
    hdul = fits.HDUList([hdu])
    out_name = f"LoTSS_{source['RA']}_{source['DEC']}.fits"
    out_path = dir + out_name
    hdul.writeto(out_path, overwrite=True)


def make_lotss_cutout(source, img, wcs):
    c = SkyCoord(source["RA"], source["DEC"], unit="deg")
    cutout = create_cutout(img, wcs, c)
    write_cutout(source, cutout)


def create_cutout(img, wcs, coord, size=5 * u.arcmin) -> Cutout2D:
    return Cutout2D(img, coord, size, wcs)


def download_ps1_from_lofar(table, outdir="./CUTOUTS/"):
    ras = table["RA"]
    decs = table["DEC"]

    ps1_table = getimages(ras, decs, filters="i")

    table["url"] = ps1_table["url"]

    download_cutouts(table)


def make_overlay(source):
    ra = source["RA"]
    dec = source["DEC"]

    lotss_cutout = f"./CUTOUTS/LoTSS_{ra}_{dec}.fits"
    ps1_cutout = f"./CUTOUTS/ps1_{ra}_{dec}.fits"

    lotss_cutout, lotss_cutout_wcs = load_image(lotss_cutout)
    ps1_cutout, ps1_cutout_wcs = load_image(ps1_cutout)

    print(lotss_cutout, lotss_cutout_wcs)

    opt_transform = ZScaleInterval()
    radio_transform = ZScaleInterval() + SqrtStretch()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ps1_cutout_wcs)

    ax.imshow(opt_transform(ps1_cutout), origin="lower")

    ax.set_autoscale_on(False)

    ax.contour(
        radio_transform(lotss_cutout),
        transform=ax.get_transform(lotss_cutout_wcs),
        levels=[1, 2, 3, 4, 5, 6],
        colors="white",
    )

    plt.show()


if __name__ == "__main__":
    radio_image = "./IMAGES/en1_radio_image.fits"

    radio_catalogue = Table.read("CATALOGUES/en1_final_component_catalogue-v1.0.fits")

    radio_catalogue.sort("Peak_flux", reverse=True)

    radio_catalogue = radio_catalogue[0:10]

    # download_ps1_from_lofar(radio_catalogue)

    en1_image, en1_wcs = load_image(radio_image)

    for source in tqdm(radio_catalogue[3:]):
        make_lotss_cutout(source, en1_image, en1_wcs)
        make_overlay(source)
