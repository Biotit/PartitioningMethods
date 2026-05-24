import numpy as np
from partitioning import Partitioning


def test_partitioning(sample_data_read):
    siteDetails = {
        "hi": 2.5,
        "zi": 4.0,
        "freq": 20,
        "length": 30,
        "PreProcessing": True,
    }
    argsOut = {"energetic_units": True, "mass_units": True, "molar_units": False}
    processing_args = {
        "density_correction": True,  # if True, density corrections are implemented during pre-processing
        "fluctuations": "LD",
        "maxGapsInterpolate": 5,  # Intervals of up to 5 missing values are filled by linear interpolation
        "RemainingData": 95,  # Only proceed with partioning if 95% of initial data is available after pre-processing
        "steadyness": False,
        "saveprocessed": False,
    }

    part = Partitioning(
        hi=siteDetails["hi"],
        zi=siteDetails["zi"],
        freq=siteDetails["freq"],
        length=siteDetails["length"],
        df=sample_data_read,
        PreProcessing=siteDetails["PreProcessing"],
        argsQC=processing_args,
        argsOut=argsOut,
    )
    part.TurbulentStats()
    part.partCEC(H=0)  # E, T, P, R
    part.partREA(H=0.0)  # E, T, P, R
    part.partCEA(H=0.0)  # E, T, P, R
    part.WaterUseEfficiency(ppath="C3")  # const_ppm, const_ratio, linear, sqrt, opt
    part.partFVS(W=part.wue["const_ppm"])  # E, T, P, R
    part.partCECw(W=part.wue["linear"])  # E, T, P, R

    output = np.array(
        [
            part.wue["const_ppm"],
            part.wue["const_ratio"],
            part.wue["linear"],
            part.wue["opt"],
            part.turbstats["LE"].magnitude,
            part.turbstats["Fc"].magnitude,
            part.fluxesCEC["Tcec"].magnitude,
            part.fluxesCEC["Pcec"].magnitude,
            part.fluxesREA["Tmrea"].magnitude,
            part.fluxesREA["Pmrea"].magnitude,
            part.fluxesCEA["Tcea"].magnitude,
            part.fluxesCEA["Pcea"].magnitude,
            part.fluxesFVS["Tfvs"].magnitude,
            part.fluxesFVS["Pfvs"].magnitude,
            part.fluxesCECw["Tcecw"].magnitude,
            part.fluxesCECw["Pcecw"].magnitude,
        ]
    )
    # With new density correction, using specific humidity for T calculation
    expected_output = np.array(
        [
            -1.05471903e-02,
            -1.19095094e-02,
            -9.26281175e-03,
            -1.95313827e-03,
            2.78027367e02,
            -9.76881156e-01,
            2.78019644e02,
            -9.76854020e-01,
            2.78027367e02,
            -9.76881156e-01,
            2.30251687e02,
            -1.21632099e00,
            2.27696304e02,
            -9.79028233e-01,
            2.67906705e02,
            -1.01164671e00,
        ]
    )
    # Old expected output with old density correction, using mixing ratio for T calculation
    # expected_output = np.array(
    #     [
    #         -1.05682305e-02,
    #         -1.19357634e-02,
    #         -9.26358302e-03,
    #         -1.89720015e-03,
    #         2.77989779e02,
    #         -9.77573590e-01,
    #         2.77989779e02,
    #         -9.77573590e-01,
    #         2.77997502e02,
    #         -9.77600746e-01,
    #         2.30218128e02,
    #         -1.21698767e00,
    #         2.27382022e02,
    #         -9.79627239e-01,
    #         2.67973775e02,
    #         -1.01198423e00,
    #     ]
    # )
    np.testing.assert_allclose(output, expected_output, rtol=1e-5, atol=1e-5)
