from astropy.io import fits
from astropy.table import Table, join_skycoord, join
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def open_table(filename, format="fits"):
    t = Table.read(filename, format=format)
    return t


def open_fits(filename):
    f = fits.open(filename)
    return f[1].data


def spec_index(flux1, flux2, freq1, freq2):
    return np.log10(flux1 / flux2) / np.log10(freq1 / freq2)


def catalogue_skycoord(t, ra, dec):
    return SkyCoord(t[ra], t[dec], unit="deg")


def cross_match_tables(c1, c2, seperation=2 * u.arcsec):
    idx, d2d, d3d = match_coordinates_sky(c1, c2)

    sep_constaint = d2d < seperation

    return sep_constaint, idx[sep_constaint]


def plot_spec_index(indeces):
    fig, ax = plt.subplots()

    ax.hist(indeces, bins=100, density=True)

    ax.vlines(-0.7, 0, ax.get_ylim()[1], color="black", linestyle="--")

    ax.set_xlabel(r"$\alpha_{\mathrm{144-3000MHz}}$")

    plt.savefig("PLOTS/spec_indeces.pdf")


def plot_2d_spec_index(index1, index2):
    fig, ax = plt.subplots()

    hist = ax.hist2d(index1, index2, bins=50, norm=LogNorm())

    plt.colorbar(hist[3], ax=ax)

    plt.plot([-5, 3], [-5, 3])

    ax.set_xlabel(r"$\alpha_{\mathrm{144-1400MHz}}$")
    ax.set_ylabel(r"$\alpha_{\mathrm{1400-3000MHz}}$")

    plt.savefig("PLOTS/spec_shape.png")


if __name__ == "__main__":
    lotss_dr3: str = "CATALOGUES/LoTSS_DR3_v0.2_srl.fits"  # Taken from https://lofar-surveys.org/dr3.html - version 0.2 - 2.7GB
    vlass: str = "CATALOGUES/VLASS_QL2_hosts.csv"  # Taken from https://cirada.ca/vlasscatalogueql0 - QL2 host ID table - 441MB
    nvss: str = "./CATALOGUES/nvss.fits"  # Taken from https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=%20NVSS - Requested unlimited as binary fits - 116MB

    try:
        matched: str = "CATALOGUES/lotss_vlass_nvss_matched.fits"
        lotss_vlass_nvss: Table = open_table(matched)
    except:
        lotss_dr3: Table = open_table(lotss_dr3)
        vlass: Table = open_table(vlass, format="csv")
        nvss: Table = open_table(nvss)

        lotss_c: SkyCoord = catalogue_skycoord(lotss_dr3, "RA", "DEC")
        vlass_c: SkyCoord = catalogue_skycoord(vlass, "RA_Host", "DEC_Host")
        nvss_c: SkyCoord = catalogue_skycoord(nvss, "RAJ2000", "DEJ2000")

        lotss_dr3["sc"] = lotss_c
        vlass["sc"] = vlass_c
        nvss["sc"] = nvss_c

        lotss_vlass: Table = join(
            lotss_dr3, vlass, join_funcs={"sc": join_skycoord(2 * u.arcsec)}
        )
        del lotss_vlass["sc_1", "sc_2"]

        lotss_c: SkyCoord = catalogue_skycoord(lotss_vlass, "RA", "DEC")
        lotss_vlass["sc"] = lotss_c

        lotss_vlass_nvss: Table = join(
            lotss_vlass, nvss, join_funcs={"sc": join_skycoord(2 * u.arcsec)}
        )
        lotss_vlass_nvss.write(
            "CATALOGUES/lotss_vlass_nvss_matched.fits", overwrite=True
        )

    a_144_1400: np.ndarray = spec_index(
        lotss_vlass_nvss["Total_flux"], lotss_vlass_nvss["S1_4"], 144e6, 1.4e9
    )
    a_1400_3000: np.ndarray = spec_index(
        lotss_vlass_nvss["S1_4"], lotss_vlass_nvss["Total_flux_source"], 1.4e9, 3e9
    )

    plot_2d_spec_index(a_144_1400, a_1400_3000)
    plot_spec_index(a_144_1400)
