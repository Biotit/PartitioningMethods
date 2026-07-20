#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from glob import glob
import multiprocessing
import datetime
from functools import partial
import logging
import os
from .Partitioning import Partitioning
from .auxfunctions import Constants, setup_logging

# Set up logger
logger = logging.getLogger(__name__)


def CallPartitioning(
    filei,
    siteDetails,
    argsQC={},
    argsOut={},
    methods=True,
    methodsWue=True,
    statistics={"TurbStats": True},
    argsQThres={},
    loadfnct="NormLoad",
    loadkwargs={},
    versatile_loadkwargs={},
):
    """
    Calls different partitioning methods and returns the data and units.

    Parameters
    ----------
    filei : str
        String with path to file being loaded and partitoned.
    siteDetails : dict
        Dictionary with details about the measurement site.

        Keys
        ----
            hi : int/float,
                Canopy mean height in meters
            zi : int/float,
                EC measurement height in meters
            freq : int/float,
                EC measurement frequency in Hz
            length : int/float,
                length of data file in minutes
            PreProcessing : bool
                If pre-processing takes place.
            ppath : str
                Type of photosynthesis ('C3' or 'C4'), for WUE calculation.
    argsQC : dict
        Contains options to be used during pre-processing regarding fluctuation extraction and if density corrections are
        necessary. All options have default values, but can be modified if needed.

        Keys
        ----
            density_correction : bool
                True if density corrections are necessary (open gas analyzer); False (closed or enclosed gas analyzer).
            fluctuations : str
                Describes the type of operation used to extract fluctuations:
                'BA': block average
                'LD': Linear detrending
                'FL': Filter low frequencies. Requires filtercut to indicate the cutoff time in minutes.
            filtercut : int
                Cutoff time in minutes for the low-pass filter. Only used if method is 'FL'.
            maxGapsInterpolate : int
                Number of consecutive gaps that will be interpolated.
            RemainingData : int
                Percentage (0-100) of the time series that should remain after pre-processing. If less than this quantity, partitioning is not implemented.
            saveprocessed : bool
                If True, the pre-processed data is saved to a CSV file in the subfolder ProcessedData.
            time_lag_correction : bool
                If True, a time lag correction is applied to the CO2 and H2O time series relative to the W time series.
            max_lag_seconds : int
                Maximum time lag in seconds to consider for correlation. Defaults to 5 seconds.
            type_lag : str
                Specifies the type of lag to consider. Options are 'positive', 'negative', or 'both'. Defaults to 'positive'.
                'Positive' means that CO2 and H2O lag behind W as expected in closed-path systems when the tube delays the signal.
            saveplotlag - bool
                If True, saves a plot of the cross-correlation function between the CO2 and H2O time series with respect to the W time series in the subfolder TimeLagCorrelationFigures.
            outfolder : str
                 If an outfolder is given the plots of the cross-correlation and the pre-processed files are saved there in the subfolders TimeLagCorrelationFigures or ProcessedData. If not, the current working directory.
    argsOut : dict
        Contains options in which units the results are given. Defaults are that the output is in mass based units.
        Possible to activate all simultanously.

        Keys
        ----
            energetic_units : bool
                True if the H2O flux shall be provided in energetic units in W/m2
            mass_units : bool
                True if the CO2 and H2O flux shall be provided as mass flux: g/(m2 s) for h2o and mg/(m2 s) for co2.
            molar_units : bool
                True if the CO2 and H2O flux shall be provided as molar flux: mmol/(m2 s)
    methods : bool or dict, default True
        If True, all available partitioning methods are used.
        If False, no method is used.
        If dict, the specified methods are used.

        Keys
        ----
        If True, the corresponding method is calculated

        MREA : bool
        CEC : bool
        CECw : bool
        CEA : bool
        FVS : bool
    methodsWue : bool or dict, default True
        If True, all available water use efficicency methods are used.
        If False, no method is used.
        If dict, the specified methods are used.

        Keys
        ----
        If True, the corresponding method is calculated

        const_ppm : bool
        const_ratio : bool
        linear : bool
        sqrt : bool
        opt : bool

    statistics : bool or dict, default {"TurbStats":True}
        If True all possible statistics are calculated.
        If False no statistics are calculated.
        If dict, the specified statistics are calculated.
        If True, the corresponding method is calculated
        Keys not used, are set to False

        Keys
        ----
        TurbStats : bool
            If True basic general turbulence statistics are calculated.
        steadyness : bool
            If True, Foken's stationarity test is implemented to check if the data is stationary.
            If False, the test is not implemented.
            The test is only informative and does not remove data, which is left to the user's discretion.
        sampledEvents : bool
            If True the time fraction and time scale of sampled events within each quadrant are calculated.

    argsQThres : dict
        Contains the quadrant thresholds stating which amount of data needs to be present within each quadrant to
        partition the fluxes.

        Keys
        ----
            cec_per_points_Q1Q2 : int
                For CEC more % of data needs to be present within quadrant 1 and 2 to partition.
                Otherwise no partitioning is performed.
            cec_per_points_each : int
                For CEC if less or at least % of data is within one of Q1 or Q2 the flux is contributed to the other quadrant.
            mrea_per_points_Q1Q2 : int
                For MREA more % of data needs to be present within quadrant 1 and 2 to partition.
                Otherwise no partitioning is performed.
            mrea_per_points_each : int
                For MREA if less or at least % of data are within one of Q1 or Q2 the flux is contributed to the other quadrant.
            cea_per_points_each : int
                For CEA more % of data needs to be in each of the necessary four quadrants Q1 and Q2 for both
                up- and downdrafts, no partitioning is performed.
            t_scale_gap_threshold : int
                For the time scale of sampled events, the minimum amount of datapoints to define a new conditionally sampled event.
            H : float or dict
                Hyperbolic threshold criteria. If not specified 0 is used for all methods.
                If float: MREA, CEC, CEA, CECw get calculated using this threshold.
                If dict:
                    Hyperbolic threshold per method used.
                    If for a method no threshold is defined its set to 0.

                    Keys
                    ----
                    MREA : float
                        Hyperbolic threshold for MREA
                    CEC : float
                        Hyperbolic threshold for CEC
                    CEA : float
                        Hyperbolic threshold for CEA
                    CECw : float
                        Hyperbolic threshold for CECw


    loadfnct : str
        Function name as a string used for loading the data.
        Available options are:
            "VersatileLoad"
                loading, renaming and recalculations can be done with this function.
                it basically makes the other loading functions useless apart from their
                shorter notation.

            "NormLoad"
                basically pd.read_csv, pass the arguments to read the data as loadkwargs.
                the index gets to be the timestamp

            "LoadBmmflux"
                custom function to read the BMMFlux high-frequency output files.
                BMMFlux is the EddyCovariance Software of the Micrometeorology Group in Bayreuth.

                See appendix of
                Thomas, C. K., Law, B. E., Irvine, J., Martin, J. G., Pettijohn, J. C., & Davis, K. J. (2009):
                    Seasonal hydrology explains interannual and seasonal variation in carbon
                    and water exchange in a semiarid mature ponderosa pine forest in central
                    Oregon. Journal of Geophysical Research: Biogeosciences, 114(G4).
                    https://doi.org/10.1029/2009JG001010

    loadkwargs : dict
        Arguments passed to the pd.read_csv() function in case of NormLoad and VersatileLoad.

    versatile_loadkwargs : dict
        Arguments passed to the VersatileLoad function.

        Keys
        ----
            timestamp_col : str or list of str, optional
                - If a list: Combines split columns (e.g., ['Year', 'Month', 'Day']) into a datetime index.
                - If a string: Converts that specific column into the datetime index.
                - If None: Converts the default DataFrame index to datetime.
            convert_gases : bool, default False
                If True, converts 'co2' and 'h2o' from mmol/m³ to mg/m³ and g/m³ respectively.
            rename_cols : dict, optional
                A dictionary mapping old column names to new ones (e.g., {"Pressure": "P"}).
            select_cols : list of str, optional
                A list of specific columns to keep. All other columns will be dropped.

    Returns
    -------
        datun : dict
            Dictionary with partitioned data and units.

        Keys
        ----
        data : dict
            Partitioned data, each key:value pair corresponds to one return value of the partioning functions.
        units : dict
            The corresponding unit, as well as key:value pair.
    """
    part_results = {}
    units = {}  # Dictionary to store units

    # Load data
    try:
        # Load which function to use for loading
        current_module = sys.modules[__name__]
        loadfnct_f = getattr(current_module, loadfnct)

        # Load data
        df = loadfnct_f(filei, loadkwargs, **versatile_loadkwargs)

    except AttributeError:
        logger.error(
            f"Error: '{loadfnct}' doesn't exist in the package! Use VersatileLoad, NormLoad or LoadBmmflux instead."
        )
        return None

    # Settings statistics
    default_stats = {"TurbStats": False, "steadyness": False, "sampledEvents": False}
    if statistics == True:
        logger.debug("Activating all statistic methods")
        # activating all methods
        statistics = {k: True for k in default_stats}
    elif statistics == False:
        logger.debug("Deactivating all statistic methods")
        statistics = default_stats
    else:
        logger.debug("Set only some statistic methods")
        statistics = {**default_stats, **statistics}
    
    # Set datetime_start for this data, BEFORE partitioning class gets
    # initialized and NaNs gets dropped
    part_results["Datetime_start"] = df.index[0]
    
    # Setting up partitioning class, including PreProcessing during init
    try:
        part = Partitioning(
            hi=siteDetails["hi"],
            zi=siteDetails["zi"],
            freq=siteDetails["freq"],
            length=siteDetails["length"],
            df=df,
            PreProcessing=siteDetails["PreProcessing"],
            argsQC=argsQC,
            argsOut=argsOut,
            sampledEventsStats=statistics["sampledEvents"],
            argsQThres=argsQThres,
        )
    except (ValueError, TypeError) as e:
        logger.error(f"Error in {filei}: {e}")
        return None

    if argsQC.get("saveprocessed"):
        # Saving pre-processed data
        metadata = {
            "Source": "PartitioningMethods",
            "RunDate": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "argsQC": str(argsQC),
            "siteDetails": str(siteDetails),
        }
        outfolder_p = argsQC.get("outfolder", "")
        path_folder = outfolder_p + "ProcessedData"
        timestamp = df.index[0].strftime("%Y%m%d-%H%M")
        path = path_folder + "/processed-%s.csv" % (timestamp)
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)
        logger.debug(f"Saving file {path}.")
        with open(path, "w") as f:
            for _setting in metadata:
                f.write(f"# {_setting}:, {metadata[_setting]}\n")
            f.write("\n")
            part.data.to_csv(f, na_rep="NaN", index=False)


    # Helper function to extract magnitude and unit
    def extract_data(source_dict, target_dict, unit_dict, suffix=""):
        for key, value in source_dict.items():
            dict_key = f"{key}{suffix}" if suffix else key
            if hasattr(value, "magnitude"):
                target_dict[dict_key] = value.magnitude
                unit_dict[dict_key] = str(value.units)  # Store the unit string
            else:
                target_dict[dict_key] = value
                unit_dict[dict_key] = ""

    # Setting hyperbolic threshold
    if "H" not in argsQThres:
        logger.debug("Set all hyperbolic thresholds to 0.")
        # if not defined, set to 0
        argsQThres["H"] = 0

    if isinstance(argsQThres["H"], (int, float)):
        logger.debug(f"Set all hyperbolic thresholds to {argsQThres['H']}.")
        # if a number set this number for all methods
        argsQThres["H"] = {
            "MREA": argsQThres["H"],
            "CEC": argsQThres["H"],
            "CECw": argsQThres["H"],
            "CEA": argsQThres["H"],
        }

    for _M in ("MREA", "CEC", "CECw", "CEA"):
        logger.debug("Some methods with hyperbolic threshold, some without.")
        if _M not in argsQThres["H"]:
            # if a method was not specified (but another was), set to 0.
            argsQThres["H"][_M] = 0

    # Calculate general turbulence statistics
    if statistics["TurbStats"]:
        logger.debug("Calculating general turbulence characteristics.")
        part.TurbulentStats()
        extract_data(part.turbstats, part_results, units)

    # Calculate steadyness statistics
    if statistics["steadyness"]:
        logger.debug("Calculating steadyness characteristics.")
        part._steadynessTest()
        extract_data(part.FokenStatTest, part_results, units)

    # Settings which method to process
    default_methods = {
        "MREA": False,
        "CEC": False,
        "CECw": False,
        "CEA": False,
        "FVS": False,
    }
    if methods == True:
        logger.debug("Activating all methods for partitioning.")
        # activating all methods
        methods = {k: True for k in default_methods}
    elif methods == False:
        logger.debug("Deactivating all methods for partitioning.")
        methods = default_methods
    else:
        logger.debug("Only some methods for partitioning activated.")
        # if a method was not named in dict set it to false
        methods = {**default_methods, **methods}

    # Processing all methods
    if methods["CEC"]:
        logger.debug("Partitioning using CEC.")
        part.partCEC(H=argsQThres["H"]["CEC"])
        extract_data(part.fluxesCEC, part_results, units)

    if methods["MREA"]:
        logger.debug("Partitioning using MREA.")
        part.partREA(H=argsQThres["H"]["MREA"])
        extract_data(part.fluxesREA, part_results, units)

    if methods["CEA"]:
        logger.debug("Partitioning using CEA.")
        part.partCEA(H=argsQThres["H"]["CEA"])
        extract_data(part.fluxesCEA, part_results, units)

    try:
        logger.debug("Calculating Water Use Efficiencies")
        # Calculate Water Use Efficicency if possible
        part.WaterUseEfficiency(ppath=siteDetails["ppath"], methodsWue=methodsWue)
        for key in part.wue.keys():
            part_results[f"W_{key}"] = part.wue[key]
            units[f"W_{key}"] = ""  # WUE usually dimensionless or handled manually
        part_results["statuswue"] = "ok"
    except ValueError as e:
        logger.error("Error while processing file %s. Error caused by: %s", filei, e)
        part_results["statuswue"] = "VPD < 0"
        for _w in part.wue.keys():
            part_results[f"statusfvs_{_w}"] = "VPD < 0. No WUE."
            part_results[f"statuscecw_{_w}"] = "VPD < 0. No WUE."
        return {"data": part_results, "units": units}

    for _w in part.wue.keys():
        # Run methods relying on water use efficency for all
        # water use efficiencies available

        if methods["FVS"]:
            logger.debug("Partitioning using FVS.")
            part.partFVS(W=part.wue[_w])
            extract_data(part.fluxesFVS, part_results, units, suffix=f"_{_w}")

        if methods["CECw"]:
            logger.debug("Partitioning using CECw.")
            part.partCECw(W=part.wue[_w], H=argsQThres["H"]["CECw"])
            extract_data(part.fluxesCECw, part_results, units, suffix=f"_{_w}")

    # Return both data and units
    datun = {"data": part_results, "units": units}
    return datun


def process(
    siteDetails,
    infolder,
    outfolder,
    loadpattern="*.csv",
    outname="PartitioningResults",
    argsQC={},
    argsOut={},
    methods=True,
    methodsWue=True,
    statistics={"TurbStats": True},
    argsQThres={},
    loadfnct="NormLoad",
    loadkwargs={},
    versatile_loadkwargs={},
    logginglevel=20,
):
    """
    Loading the raw data, (optionally) pre-process, partition and save the results.

    Parameters
    ----------
        siteDetails : dict
            Dictionary with details about the measurement site.

            Keys
            ----
                hi : int/float,
                    Canopy mean height in meters
                zi : int/float,
                    EC measurement height in meters
                freq : int/float,
                    EC measurement frequency in Hz
                length : int/float,
                    length of data file in minutes
                PreProcessing : bool
                    If pre-processing takes place.
                ppath : str
                    Type of photosynthesis ('C3' or 'C4'), for WUE calculation.

        infolder - str
            Path to the folder where the input data is located. Needs to end with slash or backslash.
        outfolder - str
            Path to the folder where the output data is located. Needs to end with slash or backslash.
        loadpattern - str, default "*.csv"
            Pattern in the filename to match for loading files.
        outname - str, default "PartitioningResults"
            Filename for the output data, excluding file ending.


        argsQC : dict
            Contains options to be used during pre-processing regarding fluctuation extraction and if density corrections are
            necessary. All options have default values, but can be modified if needed.

            Keys
            ----
                density_correction : bool
                    True if density corrections are necessary (open gas analyzer); False (closed or enclosed gas analyzer).
                fluctuations : str
                    Describes the type of operation used to extract fluctuations:
                    'BA': block average
                    'LD': Linear detrending
                    'FL': Filter low frequencies. Requires filtercut to indicate the cutoff time in minutes.
                filtercut : int
                    Cutoff time in minutes for the low-pass filter. Only used if method is 'FL'.
                maxGapsInterpolate : int
                    Number of consecutive gaps that will be interpolated.
                RemainingData : int
                    Percentage (0-100) of the time series that should remain after pre-processing. If less than this quantity, partitioning is not implemented.
                saveprocessed : bool
                    If True, the pre-processed data is saved to a CSV file in the subfolder ProcessedData.
                time_lag_correction : bool
                    If True, a time lag correction is applied to the CO2 and H2O time series relative to the W time series.
                max_lag_seconds : int
                    Maximum time lag in seconds to consider for correlation. Defaults to 5 seconds.
                saveplotlag - bool
                    If True, saves a plot of the cross-correlation function between the CO2 and H2O time series with respect to the W time series in the subfolder TimeLagCorrelationFigures.
                type_lag : str
                    Specifies the type of lag to consider. Options are 'positive', 'negative', or 'both'. Defaults to 'positive'.
                    'Positive' means that CO2 and H2O lag behind W as expected in closed-path systems when the tube delays the signal.

        argsOut : dict
            Contains options in which units the results are given. Defaults are that the output is in mass based units.
            Possible to activate all simultanously.

            Keys
            ----
                energetic_units : bool
                    True if the H2O flux shall be provided in energetic units in W/m2
                mass_units : bool
                    True if the CO2 and H2O flux shall be provided as mass flux: g/(m2 s) for h2o and mg/(m2 s) for co2.
                molar_units : bool
                    True if the CO2 and H2O flux shall be provided as molar flux: mmol/(m2 s)

        methods : bool or dict, default True
            If True, all available partitioning methods are used.
            If False, no method is used.
            If dict, the specified methods are used.

            Keys
            ----
            If True, the corresponding method is calculated

            MREA : bool
            CEC : bool
            CECw : bool
            CEA : bool
            FVS : bool

        methodsWue : bool or dict, default True
            If True, all available water use efficicency methods are used.
            If False, no method is used.
            If dict, the specified methods are used.

            Keys
            ----
            If True, the corresponding method is calculated

            const_ppm : bool
            const_ratio : bool
            linear : bool
            sqrt : bool
            opt : bool

        statistics : bool or dict, default {"TurbStats":True}
            If True all possible statistics are calculated.
            If False no statistics are calculated.
            If dict, the specified statistics are calculated.

            Keys
            ----
            If True, the corresponding method is calculated
            Keys not used, are set to False

            TurbStats : bool
                If True basic general turbulence statistics are calculated.
            steadyness : bool
                If True, Foken's stationarity test is implemented to check if the data is stationary.
                If False, the test is not implemented.
                The test is only informative and does not remove data, which is left to the user's discretion.
            sampledEvents : bool
                If True the time fraction and time scale of sampled events within each quadrant are calculated.

        argsQThres : dict
            Contains the quadrant thresholds stating which amount of data needs to be present within each quadrant to
            partition the fluxes. Also the settings regarding the hyperbolic threshold and
            about the time scale of sampled events can be given here.

            Keys
            ----
                cec_per_points_Q1Q2 : int
                    For CEC more % of data needs to be present within quadrant 1 and 2 to partition.
                    Otherwise no partitioning is performed.
                cec_per_points_each : int
                    For CEC if less or at least % of data is within one of Q1 or Q2 the flux is contributed to the other quadrant.
                mrea_per_points_Q1Q2 : int
                    For MREA more % of data needs to be present within quadrant 1 and 2 to partition.
                    Otherwise no partitioning is performed.
                mrea_per_points_each : int
                    For MREA if less or at least % of data are within one of Q1 or Q2 the flux is contributed to the other quadrant.
                cea_per_points_each : int
                    For CEA more % of data needs to be in each of the necessary four quadrants Q1 and Q2 for both
                    up- and downdrafts, no partitioning is performed.
                t_scale_gap_threshold : int
                    For the time scale of sampled events, the minimum amount of datapoints to define a new conditionally sampled event.
                H : float or dict
                    Hyperbolic threshold criteria. If not specified 0 is used for all methods.
                    If float: MREA, CEC, CEA, CECw get calculated using this threshold.
                    If dict:
                        Hyperbolic threshold per method used.
                        If for a method no threshold is defined its set to 0.

                        Keys
                        ----
                        MREA : float
                            Hyperbolic threshold for MREA
                        CEC : float
                            Hyperbolic threshold for CEC
                        CEA : float
                            Hyperbolic threshold for CEA
                        CECw : float
                            Hyperbolic threshold for CECw


        loadfnct : str
            Function name as a string used for loading the data.
            Available options are:
                "VersatileLoad"
                    loading, renaming and recalculations can be done with this function.
                    it basically makes the other loading functions useless apart from their
                    shorter notation.

                "NormLoad"
                    basically pd.read_csv, pass the arguments to read the data as loadkwargs.
                    the index gets to be the timestamp

                "LoadBmmflux"
                    custom function to read the BMMFlux high-frequency output files.
                    BMMFlux is the EddyCovariance Software of the Micrometeorology Group in Bayreuth.

                    See appendix of
                    Thomas, C. K., Law, B. E., Irvine, J., Martin, J. G., Pettijohn, J. C., & Davis, K. J. (2009):
                        Seasonal hydrology explains interannual and seasonal variation in carbon
                        and water exchange in a semiarid mature ponderosa pine forest in central
                        Oregon. Journal of Geophysical Research: Biogeosciences, 114(G4).
                        https://doi.org/10.1029/2009JG001010

        loadkwargs : dict
            Arguments passed to the pd.read_csv() function in case of NormLoad and VersatileLoad.

        versatile_loadkwargs : dict
            Arguments passed to the VersatileLoad function.

            Keys
            ----
                timestamp_col : str or list of str, optional
                    - If a list: Combines split columns (e.g., ['Year', 'Month', 'Day']) into a datetime index.
                    - If a string: Converts that specific column into the datetime index.
                    - If None: Converts the default DataFrame index to datetime.
                convert_gases : bool, default False
                    If True, converts 'co2' and 'h2o' from mmol/m³ to mg/m³ and g/m³ respectively.
                rename_cols : dict, optional
                    A dictionary mapping old column names to new ones (e.g., {"Pressure": "P"}).
                select_cols : list of str, optional
                    A list of specific columns to keep. All other columns will be dropped.
        logginglevel : int, optional, default 20
            The logging threshold for the root logger. Defaults to logging.INFO (20).
            Common values are logging.DEBUG (10), logging.INFO (20), or logging.WARNING (30).

    Saves
    ----------
       df_data : pandas.DataFrame
           Processed and partioned data as csv-file with metadata header.
    Return
    ----------
       df_data : pandas.DataFrame
           Processed and partioned data
    """
    # Start logger
    filename = f"log/run_{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}.log"
    setup_logging(outfolder + filename, level=logginglevel)
    logger.info("PartitioningMethods: Starting processing and partitioning.")
    logger.info(f"See log-file under {outfolder + filename}.")

    # Find files
    logger.info("Looking for files to process.")
    listfiles = glob(infolder + loadpattern)

    if not listfiles:
        error_msg = (
            f"No files matching pattern '{loadpattern}' were found in '{infolder}'. "
            f"Please check your path or pattern."
        )
        logger.critical(error_msg)
        raise FileNotFoundError(error_msg)

    if "outfolder" not in argsQC:
        argsQC_i = argsQC.copy()
        argsQC_i["outfolder"] = outfolder
    else:
        argsQC_i = argsQC.copy()

    # only the listfiles argument is different from run to run of the
    # CallPartitioning function, all others are held constant
    # Create a version of the function where the additional arguments are already set
    partial_CallPart = partial(
        CallPartitioning,
        siteDetails=siteDetails,
        argsQC=argsQC_i,
        argsOut=argsOut,
        methods=methods,
        methodsWue=methodsWue,
        statistics=statistics,
        argsQThres=argsQThres,
        loadfnct=loadfnct,
        loadkwargs=loadkwargs,
        versatile_loadkwargs=versatile_loadkwargs,
    )

    # Run multiprocessing
    logger.info(
        "Starting multiprocessing the files. Logging message might be mixed because of multiprocessing."
    )
    pool = multiprocessing.Pool(4)
    raw_output = pool.map(partial_CallPart, listfiles)
    pool.close()
    pool.join()

    # Filter valid results
    valid_results = [r for r in raw_output if r is not None]

    if valid_results:
        logger.debug("Found valid_results.")

        # Build the main data frame
        df_data = pd.DataFrame([r["data"] for r in valid_results])
        df_data.set_index("Datetime_start", inplace=True)
        df_data.sort_index(inplace=True)

        # Get units from the first valid result
        units_row = valid_results[0]["units"]

        # Create Metadata Header
        # We create a dictionary for the top rows
        metadata = {
            "Source": "PartitioningMethods",
            "RunDate": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "argsQC": str(argsQC),
            "siteDetails": str(siteDetails),
            "argsQThres": str(argsQThres),
            "statistics": str(statistics),
            "argsOut": str(argsOut),
        }

        # Construct Multi-index for better CSV structure
        # Mapping columns to their units
        col_units = [units_row.get(col, "") for col in df_data.columns]

        df_data = df_data.reset_index()
        column_tuples = [("", "Datetime_start")] + list(zip(col_units, df_data.columns[1:]))

        # Create a MultiIndex: Columns -> Units
        df_data.columns = pd.MultiIndex.from_tuples(column_tuples)

        # Save to CSV with Metadata prepended
        output_path = outfolder + outname + ".csv"

        with open(output_path, "w") as f:
            for _setting in metadata:
                f.write(f"# {_setting}:, {metadata[_setting]}\n")
            f.write("\n")
            df_data.to_csv(f, na_rep="NaN", index=False)

        logger.info(f"Results saved successfully to {output_path}")
        return df_data
    else:
        logger.error("No valid results. No file saved.")


# Loading functions ------------------------------------------------------------


def VersatileLoad(
    path,
    loadkwargs=None,
    timestamp_col=None,
    convert_gases=False,
    rename_cols=None,
    select_cols=None,
):
    """
    A versatile CSV loading function that accommodates standard formats,
    multi-column high-frequency timestamps, gas unit conversions, and column filtering.

    Parameters
    ----------
    path : str
        Path to the csv file to be loaded.
    loadkwargs : dict, optional
        Arguments passed directly to pd.read_csv (e.g., header, skiprows, na_values).
    timestamp_col : str or list of str, optional
        - If a list: Combines split columns (e.g., ['Year', 'Month', 'Day']) into a datetime index.
        - If a string: Converts that specific column into the datetime index.
        - If None: Converts the default DataFrame index to datetime.
    convert_gases : bool, default False
        If True, converts 'co2' and 'h2o' from mmol/m³ to mg/m³ and g/m³ respectively.
    rename_cols : dict, optional
        A dictionary mapping old column names to new ones (e.g., {"Pressure": "P"}).
    select_cols : list of str, optional
        A list of specific columns to keep. All other columns will be dropped.

    Returns
    -------
    df : pandas.DataFrame
        The loaded data.
    """

    logger.debug("Loading the data with VersatileLoad.")

    # Ensure loadkwargs is a dictionary to prevent errors if passed as None
    if loadkwargs is None:
        loadkwargs = {}

    # 1. Load the raw CSV data
    df = pd.read_csv(path, **loadkwargs)

    # 2. Process Timestamp Data
    if timestamp_col is not None:
        if isinstance(timestamp_col, list):
            # Handle multi-column split dates (like BMMFlux)
            t_present_col = [col for col in timestamp_col if col in df.columns]
            if t_present_col:
                df.index = pd.to_datetime(df[t_present_col])
        elif isinstance(timestamp_col, str):
            # Handle a single named datetime column
            if timestamp_col in df.columns:
                df.index = pd.to_datetime(df[timestamp_col])
    else:
        # Convert the default row index to datetime (like NormLoad)
        df.index = pd.to_datetime(df.index)

    # 3. Rename Columns (if provided)
    if rename_cols:
        df = df.rename(columns=rename_cols)

    # 4. Filter Specific Columns (if provided)
    if select_cols:
        # Intersect to protect against KeyError if a column is missing
        valid_cols = [col for col in select_cols if col in df.columns]
        df = df[valid_cols]

    # 5. Convert CO2 and H2O units (if requested)
    if convert_gases:
        if "co2" in df.columns:
            df["co2"] = df["co2"] * 10**-3 * Constants.MWco2.magnitude * 10**6
        if "h2o" in df.columns:
            df["h2o"] = df["h2o"] * 10**-3 * Constants.MWvapor.magnitude * 10**3

    return df


def NormLoad(path, loadkwargs, **kwargs):
    """
    Loading csv data.

    Parameters
    ----------
    path : str
        Path to the csv file to be loaded.
    loadkwargs : dict
        Arguments passed to pd.read_csv.
        Options include i.a. header, index_col, usecols, names, na_values, skiprows
    kwargs : further arguments
        ignored

    Returns
    -------
    df : pandas.DataFrame
        The loaded data.
    """
    logger.debug("Loading the data with NormLoad.")
    df = pd.read_csv(path, **loadkwargs)
    df.index = pd.to_datetime(df.index)
    return df


def LoadBmmflux(path, loadkwargs, **kwargs):
    """
    Loading csv data, optimized for BMMFlux high-frequency output files.


    Parameters
    ----------
    path : str
        Path to the csv file to be loaded.
    loadkwargs : dict
        Not used.
    **kwargs : further arguments
        ignored

    Returns
    -------
    df : pandas.DataFrame
        The loaded data.
    """
    logger.debug("Loading the data with LoadBmmflux.")
    if loadkwargs:
        logger.warning("BMMFlux loading function ignores loading arguments!")

    # Reading CSV
    df = pd.read_csv(path, header=[0], na_values=["NaN"], skiprows=[1])

    # date column adjustment and index column
    t_present_col = [
        col
        for col in ["Year", "Month", "Day", "Hour", "Minute", "Second", "Millisecond"]
        if col in df.columns
    ]
    df.index = pd.to_datetime(df[t_present_col])

    # Rename columns and drop unused
    nec_cols = ["u", "v", "w", "Ts", "co2", "h2o", "P"]
    df = df.rename(columns={"Pressure": "P"})[nec_cols]

    # Recalculations
    # u, v, w, Ts and P already in correct units
    # need to recalculate co2 and h2o from mmol/m3 in mg/m3 and g/m3
    df["co2"] = df["co2"] * 10**-3 * Constants.MWco2.magnitude * 10**6
    df["h2o"] = df["h2o"] * 10**-3 * Constants.MWvapor.magnitude * 10**3

    return df
