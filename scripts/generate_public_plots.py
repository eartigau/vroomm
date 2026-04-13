from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def snr_curve_vs_magnitude(
    r_mag: np.ndarray,
    exptime: float,
    resolution: float = 120000,
    r_primary_cm: float = 160,
    wave0_angstrom: float = 6204.29,
    efficiency: float = 0.1,
    read_noise: float = 3.0,
    dark_current: float = 1.6e-4,
    spatial_extent_pix: float = 20,
) -> np.ndarray:
    """Estimate SNR per resolution element for a simple photon-noise model."""
    area = np.pi * (r_primary_cm / 2.0) ** 2
    zp = 2.49767e-9 * area
    res_element = wave0_angstrom / resolution
    flux_erg = zp * res_element * exptime * 10 ** (-r_mag / 2.5) * efficiency
    photon_energy_erg = 1.602176634e-12 * (12398.418 / wave0_angstrom)
    nphot = flux_erg / photon_energy_erg
    noise = np.sqrt(nphot + read_noise**2 * spatial_extent_pix + dark_current * spatial_extent_pix * exptime)
    return nphot / noise


def simple_photon_snr(r_mag: np.ndarray, t_total: float, exptime: float, resolution: float, efficiency: float) -> np.ndarray:
    """Approximate integrated SNR scaling adapted from simu_photo_rate.py assumptions."""
    r_primary_cm = 160
    wave0_angstrom = 6204.29
    spatial_extent = 10
    spectral_extent = 3

    area = np.pi * (r_primary_cm / 2.0) ** 2
    zp = 2.49767e-9 * area
    res_element = wave0_angstrom / resolution
    flux_erg = zp * res_element * exptime * 10 ** (-r_mag / 2.5) * efficiency
    photon_energy_erg = 1.602176634e-12 * (12398.418 / wave0_angstrom)
    nphot = flux_erg / photon_energy_erg
    nphot_per_pix = nphot / spatial_extent

    # Total SNR from short exposures accumulated over t_total.
    return np.sqrt((nphot_per_pix / exptime) * t_total * spatial_extent * 3 * spectral_extent)


def make_plots(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    r_mag = np.arange(5, 26, 1)

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    colors = {300: "#1f77b4", 900: "#2ca02c", 1800: "#d62728"}

    for exptime in [300, 900, 1800]:
        snr = snr_curve_vs_magnitude(r_mag, exptime=exptime)
        ax.plot(r_mag, snr, label=f"{exptime}s", color=colors[exptime], linewidth=2)

    ax.set_yscale("log")
    ax.set_xlim(5, 25)
    ax.set_ylim(0.2, 700)
    ax.axhline(200, color="k", linestyle="--", alpha=0.6, label="~1 m/s regime")
    ax.axhline(40, color="k", linestyle="-.", alpha=0.6, label="~5 m/s regime")
    ax.axhline(1, color="k", linestyle=":", alpha=0.6, label="~100 m/s regime")
    ax.set_xlabel("r magnitude")
    ax.set_ylabel("SNR per resolution element")
    ax.set_title("VROOMM SNR scaling (analytic)")
    ax.grid(alpha=0.2)
    ax.legend(ncol=2)
    fig.savefig(outdir / "vroomm_snr_public.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    r_mag_dense = np.linspace(8, 21, 120)
    configs = [
        (600, 0.10, 120000, "10 min, eff=0.10, R=120k", "#1f77b4"),
        (1800, 0.10, 120000, "30 min, eff=0.10, R=120k", "#ff7f0e"),
        (1800, 0.15, 80000, "30 min, eff=0.15, R=80k", "#2ca02c"),
    ]

    for t_total, eff, res, label, color in configs:
        snr = simple_photon_snr(r_mag_dense, t_total=t_total, exptime=0.25, resolution=res, efficiency=eff)
        ax.plot(r_mag_dense, snr, label=label, color=color, linewidth=2)

    ax.set_yscale("log")
    ax.set_xlim(8, 21)
    ax.set_ylim(2, 2e4)
    ax.set_xlabel("r magnitude")
    ax.set_ylabel("Integrated SNR")
    ax.set_title("Photon-rate driven SNR scenarios")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.savefig(outdir / "vroomm_photon_snr_grid.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    make_plots(root / "figures")
    print("Saved plots to figures/")
