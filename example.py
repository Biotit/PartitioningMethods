#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Explaining the possibilities of the package PartitioningMethods
"""

# Load package ----------------------------------------------------------------
# Set working directory to location of this file
# then possible to load module without installing
import src.partitioning as partitioning

# If installed (see README.md) just use
# import partitioning


# Run example dataset and partition using all methods -------------------------

# Running like this producing the file PartitioningResults_example.csv
data = partitioning.process(
    # ----
    # Mandatory settings
    # ----
    # 1) Information about the measurement site and data including units:
    # needs to be modified to your case
    siteDetails={
        "hi": 2.5,  # Canopy mean height in meters
        "zi": 4.0,  # EC measurement height in meters
        "freq": 20,  # EC measurement frequency in Hz
        "length": 30,  # length of data file in minutes
        "PreProcessing": True,  # If True, input is raw data and pre-processing
        # is applied before partitioning (e.g. density corrections
        # and fluctuations are computed)
        "ppath": "C3",  # Main plant photosynthesis type, for calculation of WUE.
        # Options: "C3" or "C4".
    },
    # 2) path to data input folder, needs to END WITH SLASH!
    infolder="RawData30min/",
    # 3) path to the folder where the results are saved, needs to END WITH SLASH!
    outfolder="results/",
    # ----
    # Optional settings
    # ----
    # 4) Pattern in the filename to match for loading files. Default: "*.csv"
    loadpattern="*.csv",
    # 5) Output name
    outname="PartitioningResults_example",
    # 6) Information about the data processing
    # ----
    # PLEASE MAKE SURE YOU HAVE THE PROPER SETTINGS HERE
    # ----
    # All options have default values assuming that raw data is provided and measured by an open-path gas analyzer and a 3D sonic anemometer
    # For usual usecase better deactivate saveprocessed to save data and time.
    # Also think about which corrections you need.
    # For this example we activate everything to get every output.
    argsQC={
        "physical_bounds": True,  # If True, data outside of specified physical bounds is set to NA.
        "despike": True,  # If True, outliers in the data are removed by despiking.
        "coord_rotation": True,  # If True, a double coordinate rotation is performed to deminish the mean vertical wind speed.
        "density_correction": True,  # If True, density corrections are implemented during pre-processing (depends on type of gas analyzer used). For closed-path gas analyzers, set "density_correction" to False
        "fluctuations": "LD",  # If "LD", linear detrending is applied to the data. BA (block averaging) and FL (filter low freqencies) are also available
        "filtercut": 5,  # Cutoff timescale to filter low frequencies (in minutes). Needed when FL is selected as fluctuation method
        "maxGapsInterpolate": 5,  # Intervals of up to 5 missing values are filled by linear interpolation
        "RemainingData": 95,  # Only proceed with partioning if 95% of initial data is available after pre-processing
        "saveprocessed": True,  # If True, save the intermediate processed data including all corrections and fluctuations in the outfolder
        "time_lag_correction": True,  # If True, a time lag correction is applied to the CO2 and H2O time series relative to the W time series
        "max_lag_seconds": 5,  # Maximum time lag in seconds to consider for cross correlation analyses
        "saveplotlag": True,  # If True, saves a plot of the cross-correlation function between the CO2 and H2O time series with respect to the W time series in the outfolder
        "type_lag": "both",  # Specifies the type of lag to consider ('negative', 'positive', or 'both')
        "UnitBorders": {  # define data range in between the median of the data has to be, otherwise an error is raised
            "Ts": (0, 70),  # C
            "co2": (200, 1500),  # mg/m3
            "h2o": (0, 50),  # g/m3
            "P": (60, 150),  # kPa
        },
        "PhysicalBounds": {  # define data range in between the values have to be, otherwise the individual values are set to NaN
            "u": (-20, 20),  # m/s
            "v": (-20, 20),  # m/s
            "w": (-20, 20),  # m/s
            "Ts": (-10, 50),  # Celsius
            "co2": (0, 1500),  # mg/m3
            "h2o": (0, 40),  #  g/m3
            "P": (60, 150),  # kPa
        },
    },
    # 7) Units in output. Defaults are that the output is in mass based units.
    argsOut={
        "energetic_units": True,  # return latent heat flux in energetic units (W/m2)
        "mass_units": True,  # return water and carbon fluxes in mass units: g/(m2 s) for h2o and mg/(m2 s) for co2
        "molar_units": True,  # return water and carbon fluxes in molar units: both in mmol/(m2 s)
    },
    # 8) Which methods to run. Default all activated
    methods={"MREA": True, "CEC": True, "CECw": True, "CEA": True, "FVS": True},
    # 9) Which Methods to calculate leaf-level
    # water use efficicency for CECw and FVS. Default all activated.
    # For C4 plants sqrt and opt are not available.
    methodsWue={
        "const_ppm": True,
        "const_ratio": True,
        "linear": True,
        "sqrt": True,
        "opt": True,
    },
    # 10) Which statistics to calculate?
    # If True, all get calculated.
    # Default default {"TurbStats":True}.
    statistics={
        "TurbStats": True,  # General Turbulence statistics
        "steadyness": True,  # Steadyness Test
        "sampledEvents": True,  # Statistics about time fraction and time scale of sampled events
    },
    # 11) Settings for the quadrant thresholds, time scale of sampled events and
    # hyperbolic thresholds.
    # If not sure what to do, just dont specify (empty dictionary {})
    # and let the defaults do the rest.
    argsQThres={
        "cec_per_points_Q1Q2": 15,  # CEC: smallest percentage of points that must be available in the first two quadrants
        "cec_per_points_each": 3,  # CEC: smallest percentage of points in each quadrant
        "cecw_per_points_Q1Q2": 0,  # CECw: smallest percentage of points that must be available in the first two quadrants
        "cecw_per_points_each": 0,  # CECw: smallest percentage of points in each quadrant
        "mrea_per_points_Q1Q2": 15,  # MREA: smallest percentage of points that must be available in the first two quadrants
        "mrea_per_points_each": 3,  # MREA: smallest percentage of points in each quadrant
        "cea_per_points_Q1Q2": 0,  # CEA: smallest percentage of points that must be available in the first two quadrants in both up and downdrafts
        "cea_per_points_each": 0,  # CEA: smallest percentage of points in each quadrant in both up and downdrafts
        "t_scale_gap_threshold": 10,  # For the time scale of sampled events, the minimum amount of datapoints to define a new conditionally sampled event (Thomas et al. 2008)
        "H": {  # Hyperbolic threshold criteria. If not specified 0 is used for all methods.
            "MREA": 0.25,  # only threshold for MREA, and its 0.25, can define also for other methods
        },  # Otherwise if no dict but float: MREA, CEC, CEA, CECw get calculated using this threshold.
    },
    # ----
    # Loading data:
    #
    # Ensure that the header of the file constains the following variables:
    # "date","u","v","w","Ts","co2","h2o","Tair","P"
    # Before processing the data,
    # ensure that the units are correct (convert if necessary)
    # "date": [yyyy-mm-dd HH:MM:SS]
    #  "u","v","w": [m/s]
    #         "Ts": [oC]
    #        "co2": [mg_CO2/m3]
    #        "h2o": [g_h2o/m3]
    #       "Tair": [oC]
    #          "P": [kPa]
    # If you have different input data, adjust using
    # VersatileLoad and versatile_loadkwargs
    # see their settings for recalculations of the units,
    # combination of several columns to create a proper timestamp
    # or renaming if your column names dont match to those listed above.
    # ----
    # 12) Which function to load the data?
    # VersatileLoad is the most advanced, however, for just loading basic files
    # without recalculations etc.
    # NormLoad can be good as well.
    loadfnct="NormLoad",
    # 13) Options for loading the data directly passed to the
    # default pandas function: pd.read_csv()
    loadkwargs={
        "header": None,
        "index_col": 0,
        "usecols": [0, 1, 2, 3, 4, 6, 8, 11, 12],
        "names": ["date", "u", "v", "w", "Ts", "co2", "h2o", "Tair", "P"],
        "na_values": ["NAN", -9999, "-9999", "#NA", "NULL"],
        "skiprows": [0],
    },
    # 14) Further loading options for VersatileLoad
    # in our case we dont need any. See documentation for available options.
    versatile_loadkwargs={},
    # 15) Logger level
    logginglevel=20,  # corresponds to logger.INFO, which is also the default
    # to see all logging for development go to 10 (logger.DEBUG).
    # A logger file is created in the log- folder in the outfolder path.
)

# Just running as script -----------------------------------------------------

# If we dont want to work with an IDE etc., we can just run the
# script from the terminal using
# python example.py

# BUT: We need to have all required packages installed in the environment.
# For further information see the README.md


# For experts: Go in depth into the package ----------------------------------

# We can also only run one dataset and one method with the function process()
# Just specify only one method and only one dataset in the input folder.
# Thats the easiest way!

# But if you want to understand how the code works, we can also do this:

# 1) Load the data.
# We can do that with the functions from the package or on our own.
# But first we need to know which file to run:

# we can use glob to see all files in the folder that are csv.
from glob import glob

listfiles = glob("RawData30min/*.csv")

# then lets load the data, we can use the functions from the package
# or our own.
loadkwargs = {
    "header": None,
    "index_col": 0,
    "usecols": [0, 1, 2, 3, 4, 6, 8, 11, 12],
    "names": ["date", "u", "v", "w", "Ts", "co2", "h2o", "Tair", "P"],
    "na_values": ["NAN", -9999, "-9999", "#NA", "NULL"],
    "skiprows": [0],
}
df = partitioning.wrapper.NormLoad(path=listfiles[0], loadkwargs=loadkwargs)
df

# 2) Lets create the Partitioning object class.
# First define the siteDetails as before.
siteDetails = {
    "hi": 2.5,  # Canopy mean height in meters
    "zi": 4.0,  # EC measurement height in meters
    "freq": 20,  # EC measurement frequency in Hz
    "length": 30,  # length of data file in minutes
    "PreProcessing": True,  # If True, input is raw data and pre-processing
    # is applied before partitioning (e.g. density corrections
    # and fluctuations are computed)
    "ppath": "C3",  # Main plant photosynthesis type, for calculation of WUE.
    # Options: "C3" or "C4".
}

# Then creating the object of class partitioning
part = partitioning.Partitioning(
    # Please look at the documentation to find what all these arguments do.
    hi=siteDetails["hi"],
    zi=siteDetails["zi"],
    freq=siteDetails["freq"],
    length=siteDetails["length"],
    df=df,
    PreProcessing=siteDetails["PreProcessing"],
    # For all further arguments we just use the defaults here
    # But PLEASE DEFINE THEM PROPERLY WHEN RUNNING SERIOUSLY.
)

# 3) We can now use the methods of this class to partition the data.
part.TurbulentStats()
# And we can access the data using the respective dictionaries.
part.turbstats

# Same for accessing the methods
part.partCEC()
part.fluxesCEC
# with .magnitude we can get rid of the unit
part.fluxesCEC["Ecec_m"].magnitude

# Save the results in a dictionary (use .magnitude to get the value and .units to get the units)
part_results = {}
units = {}
for key, value in part.fluxesCEC.items():
    dict_key = f"{key}"
    if hasattr(value, "magnitude"):
        part_results[dict_key] = value.magnitude
        units[dict_key] = str(value.units)  # Store the unit string
    else:
        part_results[dict_key] = value
        units[dict_key] = ""
# Now we have the values in this dictionary
part_results
# And the units in this
units

# See the documentation of the class Partitioning for further details.
# Also ... at this point just look the source code :)
