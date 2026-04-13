"""Explicit analytic SNR-vs-magnitude model for baseline VROOMM scenarios.

This script is intentionally verbose and unit-explicit for new team members.
It computes photons per resolution element and then propagates detector noise.
"""

import matplotlib.pyplot as plt
import numpy as np


# -------------------------------
# User-facing configuration block
# -------------------------------

# Magnitude grid used for SNR curves.
R_MAG_GRID = np.arange(5, 28)

# Exposure times to compare (seconds).
EXPOSURE_TIMES_S = [1800]

# Telescope/instrument assumptions.
PRIMARY_DIAMETER_CM = 160
CENTRAL_WAVELENGTH_A = 6204.29
RESOLUTION = 120000
THROUGHPUT = 0.1

# Detector assumptions.
SPATIAL_EXTENT_PIX = 20
READ_NOISE_E = 3.0
DARK_CURRENT_OMM_E_PER_S_PER_PIX = 1.6e-4
DARK_CURRENT_4M_E_PER_S_PER_PIX = 2.9 / 3600.0

# Optional guiding diagnostic settings.
GUIDING_BANDPASS_FRACTION = 4
GUIDING_THROUGHPUT = 0.04


def photons_per_resolution_element(r_mag: np.ndarray, exptime_s: float) -> np.ndarray:
    """Return photon counts per resolution element for each magnitude in r_mag.

    Uses a zero-point flux in erg/cm^2/s/A and converts to photon counts at
    CENTRAL_WAVELENGTH_A.
    """
    collecting_area_cm2 = np.pi * (PRIMARY_DIAMETER_CM / 2.0) ** 2
    zero_point_erg_per_cm2_s_a = 2.49767e-9

    # Flux entering telescope per angstrom.
    zero_point_erg_per_s_a = zero_point_erg_per_cm2_s_a * collecting_area_cm2

    # Width of one resolution element in angstrom.
    delta_lambda_a = CENTRAL_WAVELENGTH_A / RESOLUTION

    # Total energy per resolution element collected during exptime_s.
    signal_erg = (
        zero_point_erg_per_s_a
        * delta_lambda_a
        * exptime_s
        * 10 ** (-r_mag / 2.5)
        * THROUGHPUT
    )

    # Photon energy at CENTRAL_WAVELENGTH_A in erg.
    photon_energy_erg = 1.602176634e-12 * (12398.418 / CENTRAL_WAVELENGTH_A)
    return signal_erg / photon_energy_erg


def snr_with_read_noise(n_phot: np.ndarray, exptime_s: float, dark_current: float) -> np.ndarray:
    """SNR including photon noise, read noise, and dark current."""
    variance = n_phot + READ_NOISE_E**2 * SPATIAL_EXTENT_PIX + dark_current * SPATIAL_EXTENT_PIX * exptime_s
    return n_phot / np.sqrt(variance)


def snr_emccd_like(n_phot: np.ndarray, exptime_s: float, dark_current: float) -> np.ndarray:
    """SNR without read-noise term (EMCCD-like approximation)."""
    variance = n_phot + dark_current * SPATIAL_EXTENT_PIX * exptime_s
    return n_phot / np.sqrt(variance)


def main() -> None:
    colors = ["red", "orange", "green", "blue"]
    fig, ax = plt.subplots(layout="constrained")

    for i, exptime_s in enumerate(EXPOSURE_TIMES_S):
        n_phot = photons_per_resolution_element(R_MAG_GRID, exptime_s)

        # Optional guiding diagnostic retained from original script.
        if exptime_s == 60:
            n_phot_per_sec = n_phot / exptime_s
            guiding_photons = (
                n_phot_per_sec
                * (RESOLUTION / GUIDING_BANDPASS_FRACTION)
                * GUIDING_THROUGHPUT
                / THROUGHPUT
            )
            for mag, gphot in zip(R_MAG_GRID, guiding_photons):
                print(f"guiding photons per second {gphot:.0f}, mag = {mag:.0f}")

        snr_omm = snr_with_read_noise(n_phot, exptime_s, DARK_CURRENT_OMM_E_PER_S_PER_PIX)
        snr_omm_emccd = snr_emccd_like(n_phot, exptime_s, DARK_CURRENT_OMM_E_PER_S_PER_PIX)

        # Keep original 4-m comparison scaling (x6 photons) for continuity.
        n_phot_4m = n_phot * 6
        snr_4m = snr_with_read_noise(n_phot_4m, exptime_s, DARK_CURRENT_4M_E_PER_S_PER_PIX)

        color = colors[i % len(colors)]
        ax.plot(R_MAG_GRID, snr_omm, label=f"{exptime_s} s - OMM", color=color)
        ax.plot(R_MAG_GRID, snr_omm_emccd, label=f"{exptime_s} s - OMM, EMCCD", color=color, linestyle="--")
        ax.plot(R_MAG_GRID, snr_4m, label=f"{exptime_s} s - 4-m telescope", color=color, linestyle="dashdot")

    ax.set_xlabel("r magnitude")
    ax.set_ylabel("SNR per resolution element")
    ax.set_yscale("log")

    # RV precision guide rails.
    ax.plot([5, 30], [1, 1], "k:", alpha=0.5, label="100 m/s")
    ax.plot([5, 30], [200, 200], "k--", alpha=0.5, label="1 m/s")
    ax.plot([5, 30], [40, 40], "k-.", alpha=0.5, label="5 m/s")

    ax.set(ylim=[0.1, 600], xlim=[5, 26], title="VROOMM SNR")
    ax.legend()
    plt.savefig("vroomm_snr.png")
    plt.show()


if __name__ == "__main__":
    main()
