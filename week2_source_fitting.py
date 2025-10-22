from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from tqdm import tqdm
from astropy.visualization import ZScaleInterval, PercentileInterval, SqrtStretch
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from photutils.background import (
    Background2D,
    MedianBackground,
    LocalBackground,
    MMMBackground,
)
from astropy.convolution import convolve
from photutils.segmentation import (
    detect_sources,
    deblend_sources,
    SourceFinder,
    SourceCatalog,
    make_2dgaussian_kernel,
)
from photutils.psf import IterativePSFPhotometry, SourceGrouper, GaussianPSF
import numpy as np


def load_fits(path):
    with fits.open(path) as f:
        data = f[0].data.squeeze().squeeze()
        wcs = WCS(f[0].header).celestial
        return data, wcs


def load_table(path, format="fits"):
    return Table.read(path, format=format)


def sort_table(t, column="Total_flux"):
    t.sort(column, reverse=True)


def load_cutout(source):
    ra = source["RA"]
    dec = source["DEC"]
    path = f"./CUTOUTS/LoTSS_{ra}_{dec}.fits"
    data, wcs = load_fits(path)
    return data, wcs


def ax_show_image(data, ax):
    transform = ZScaleInterval() + SqrtStretch()

    trans_data = transform(data)

    ax.imshow(trans_data, aspect="equal", origin="lower")


def fit_image(data):
    return data


def fit_source(source, show_plot=False):
    data, wcs = load_cutout(source)

    bkg_estimator = MedianBackground()

    bkg = Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)

    data -= bkg.background

    threshold = 5.0 * bkg.background_rms

    kernel = make_2dgaussian_kernel(3.0, size=5)

    convolved_data = convolve(data, kernel)

    finder = SourceFinder(npixels=10, progress_bar=False)
    segment_map = finder(convolved_data, threshold)

    cat = SourceCatalog(data, segment_map, convolved_data=convolved_data)
    # print(cat)

    columns = [
        "label",
        "xcentroid",
        "ycentroid",
        "area",
        "segment_flux",
        "kron_flux",
        "max_value",
    ]
    tbl = cat.to_table(columns=columns)
    # print(tbl)

    if show_plot:
        fig, axs = plt.subplots(1, 3)

        ax_show_image(data, axs[0])

        axs[1].imshow(segment_map, origin="lower", cmap=segment_map.cmap)

        cat.plot_kron_apertures(ax=axs[0], color="white", lw=1.5)
        cat.plot_kron_apertures(ax=axs[1], color="white", lw=1.5)

        plt.show()

    return tbl


def fit_psf_source(source, show_plot=False):
    data, wcs = load_cutout(source)

    coord = SkyCoord(source["RA"], source["DEC"], unit="deg")

    x, y = wcs.world_to_pixel(coord)


def sum_fluxes(tbl, n_sources=3):
    # Needs a beam correction down the line
    tbl.sort("kron_flux", reverse=True)

    sub_tbl = tbl[0:n_sources]

    return sum(sub_tbl["kron_flux"])


def get_peak_flux(tbl):
    return max(tbl["max_value"])


def get_pixels_per_beam(image):
    h = fits.getheader(image)

    bmaj = h.get("BMAJ")
    bmin = h.get("BMIN")

    bmaj_arcsec = float(bmaj) * 3600.0
    bmin_arcsec = float(bmin) * 3600.0

    px1 = abs(h["CDELT1"]) * 3600.0
    px2 = abs(h["CDELT2"]) * 3600.0
    pixscale = np.sqrt(px1 * px2)

    A_beam = (np.pi / (4.0 * np.log(2.0))) * bmaj_arcsec * bmin_arcsec  #
    A_pix = pixscale**2
    Npix_per_beam = A_beam / A_pix

    return Npix_per_beam


def compare_flux_plot(fitted_fluxes, peak_fluxes):
    fig, axs = plt.subplots(2, 2)

    axs = axs.ravel()

    axs[0].set_title("Total Fluxes")
    axs[0].plot([1e-1, 1e1], [1e-1, 1e1], color="k", linestyle="--")
    axs[0].scatter(edfn[0:200]["Total_flux"], fitted_fluxes / n_pix_per_beam, s=2.0)

    axs[0].set_xlabel("Real Flux [Jy]")
    axs[0].set_ylabel("Fitted Flux [Jy]")

    axs[0].set_xscale("log")
    axs[0].set_yscale("log")

    axs[1].set_title("Peak Fluxes")

    axs[1].plot([1e-2, 1e1], [1e-2, 1e1], color="k", linestyle="--")

    axs[1].scatter(edfn[0:200]["Peak_flux"], peak_fluxes, s=2.0)

    axs[1].set_xlabel("Real Flux [Jy]")
    axs[1].set_ylabel("Fitted Flux [Jy]")

    axs[1].set_xscale("log")
    axs[1].set_yscale("log")

    delta_total = 100 * (edfn[0:200]["Total_flux"] - (fitted_fluxes / n_pix_per_beam)) / edfn[
        0:200
    ]["Total_flux"] 
    axs[2].scatter(edfn[0:200]["Total_flux"], delta_total, s= 2.0)

    axs[2].set_xscale("log")
    axs[2].set_yscale('symlog')

    delta_peak = 100 * (edfn[0:200]["Peak_flux"] - peak_fluxes) / edfn[0:200]["Peak_flux"]

    axs[3].scatter(edfn[0:200]["Peak_flux"], delta_peak, s = 2.0)

    axs[3].set_xscale("log")
    axs[3].set_yscale("symlog")

    plt.tight_layout()

    plt.savefig("./PLOTS/flux_measurements.pdf")


if __name__ == "__main__":
    edfn = load_table("CATALOGUES/LOFAR_EDFN_v1.0.fits")

    sort_table(edfn)

    fitted_fluxes = []
    peak_fluxes = []

    n_pix_per_beam = get_pixels_per_beam(
        "IMAGES/image_full_ampphase_di_m.NS.int.restoredcorrectedflux.fits"
    )

    for source in tqdm(edfn[0:200]):
        tbl = fit_source(source, show_plot=False)

        total_flux = sum_fluxes(tbl, n_sources=3)
        fitted_fluxes.append(total_flux)
        peak_fluxes.append(get_peak_flux(tbl))

    compare_flux_plot(fitted_fluxes, peak_fluxes)
