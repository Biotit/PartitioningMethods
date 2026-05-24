
<!-- hidden DOI: [![DOI](https://zenodo.org/badge/441544177.svg)](https://zenodo.org/badge/latestdoi/441544177) -->

---
# Contact information 

Original author: Einara Zahn\
email: einaraz@princeton.edu, einara.zahn@gmail.com

This version: Daniel Schöndorf
email: Daniel.Schoendorf@uni-bayreuth.de

## Flux Partitioning Package

This Python package implements five partitioning methods to separate evapotranspiration (ET) and CO<sub>2</sub> fluxes (Fc) into ground (evaporation and respiration) and plant (transpiration and net CO<sub>2</sub> assimilation) fluxes. It processes instantaneous raw eddy covariance measurements and returns fluxes, flux components, and additional turbulent variables.

The package includes pre-processing steps such as coordinate rotation and density corrections. While it does not perform all post-processing corrections provided by other eddy covariance flux software, users are encouraged to reach out with inquiries about adding additional pre- and post-processing techniques or to contribute to the code.

## New in modified version (3.0.0)

The modified version includes
- each pre-processing step can be activated seperately.
- the time fraction and time scale of sampled events by each partitioning method can be calculated.
- the sampling thresholds per quadrant can be set manually for each method.
- the output can be set to mass based units, molar units or energetic units.
- if timestamps are missing in the input data, they are filled with NaN to ensure a continous dataset.
- wrapper functions automatically performing the partitioning and save the results in a csv-file with metadata header and unit information for each quantity. Detailed settings enable to run for almost any use case.
- further data loading functions are included.
- logger files for each run are saved in the folder log.
- minor bug fixes in density correction and C4 plant water use efficiency calculation.

<!--  See the [Documentation](https://einaraz.github.io/PartitioningMethods/) here. -->
See the [Documentation](https://biotit.github.io/PartitioningMethods/) here.


# Installation Guide

To install the package, follow these steps:

1. **Clone or download the repository** to your computer.

    ```sh
    git clone https://github.com/Biotit/PartitioningMethods
    ```

2. **Navigate to the main folder** of the repository:

    ```sh
    cd PartitioningMethods/
    ```

3. **Create virtual environments** (not required, but recommended to avoid conflicts with local python packages)

**A) Create and activate a virtual environment using base python**:

    > **Note**: A virtual environment is an isolated environment in which you can install packages without affecting your system-wide Python installation. This helps prevent conflicts between package versions and keeps your project dependencies organized. Alternatively, use conda for environment management. See step 3.B) below.

    - **Create the virtual environment:**

      This command creates a directory named `Partitioning` that contains the virtual environment. You only need to do this once per project.

      ```sh
      python -m venv Partitioning
      ```

    - **Activate the virtual environment:**

      Activation adjusts your environment to use the packages installed in the virtual environment instead of the global Python installation. You need to activate the environment every time you start a new terminal session or work on your project.

      - **On macOS/Linux:**

        ```sh
        source Partitioning/bin/activate
        ```

      - **On Windows:**

        ```sh
        Partitioning\Scripts\activate
        ```

      After activation, your command prompt should change to show the name of the virtual environment, typically `(Partitioning)`.

      If you see `(Partitioning)` at the beginning of your command prompt, the virtual environment is active. If you don't see it, try running the activation command again. 

**B) Create and activate a virtual environment using conda**:

    > **Note**: Conda is a package manager with several advantages over pip. If you use conda already, it can be useful here.

    - **Create the virtual environment:**

      This command creates a directory named `Partitioning` that contains the virtual environment. You only need to do this once per project. In this environment pip is installed for easy installation of the package.

      ```sh
      conda create -n Partitioning pip
      ```

    - **Activate the virtual environment:**

      Activation adjusts your environment to use the packages installed in the virtual environment instead of the global Python installation. You need to activate the environment every time you start a new terminal session or work on your project.
      
        ```sh
        conda activate Partitioning
        ```
      After activation, your command prompt should change to show the name of the virtual environment, typically `(Partitioning)`.

      If you see `(Partitioning)` at the beginning of your command prompt, the virtual environment is active. If you don't see it, try running the activation command again. 


4. **Install the package** using pip:

    With the virtual environment activated, install your package:

    ```sh
    pip install .
    ```

    This command installs the package partitioning and all other python dependencies in the virtual environment.
   
6. **Verify the installation** by importing the package in a Python interpreter:

    To check if the package is installed correctly, open a Python interpreter and try importing it:

    ```sh
    python -c "import partitioning"
    ```

    If no errors occur, the package is installed correctly.


### Additional Information

- **Deactivating the Virtual Environment:**

  When you’re done working in the virtual environment, you can deactivate it by running:

  ```sh
  deactivate
  ```
  When using conda use
  ```sh
  conda deactivate
  ```
  
  **Important:** After installing the package `partitioning` inside the virtual environment, you must activate the `Partitioning` environment each time you work on the project. This ensures that you are using the correct package versions and dependencies specified for your project.


## Format of Input Text Files

This script works with text files separated by commas (CSV format). It reads eddy-covariance time series of any length, although 30-minute intervals are typically used under neutral to unstable conditions.

The following variables are required by the code when using raw high-frequency data as input (see the code for requirements when using pre-processed data as input):

- **index**: Date and time of acquisition in the format `[yyyy-mm-dd HH:MM]`.
- **w**: Vertical velocity component (m/s).
- **u**: Streamwise velocity component (m/s).
- **v**: Cross-stream velocity component (m/s).
- **Ts**: Sonic temperature (Celsius).
- **co2**: Carbon dioxide density (mg_CO2/m³).
- **h2o**: Water vapor density (g_H2O/m³).
- **P**: Pressure (kPa).

If your columns have different names or units or the datetime information is formatted differently, you probably can use the option `VersatileLoad` and `versatile_loadkwargs` in the process()- function to adjust to your data. See the code documentation for this.
  
## How to use it
For a complete example of how to use the partitioning module, see ```example.py```.
<!-- See the [Documentation](https://einaraz.github.io/PartitioningMethods/) for more details. -->
See the [Documentation](https://biotit.github.io/PartitioningMethods/) for more details.

The easiest is the usage of the `process` function:
```python
  import partitioning
  
  data = partitioning.process(
  	 # specifiy your custom arguments here! See example.py file!
  )
  ```

When using the function `process`, in the `outfolder` an output file with the results is created and a folder `log` with log files for all runs.

---
## Description

The package provides tools for processing and analyzing high-frequency eddy-covariance data. Key features include:

1. **Quality Control and Pre-processing**:
   - Checks for physically realistic values
   - Detects outliers
   - Implements time lag corrections
   - Rotates coordinates
   - Extracts fluctuations (block average, linear detrending, and filtering operations available)
   - Applies density corrections for CO<sub>2</sub> and H<sub>2</sub>O measured by open-path gas analyzers ("instantaneous" WPL correction, based on the paper [Detto and Katul, 2007](https://link.springer.com/article/10.1007%2Fs10546-006-9105-1))
   - Performs stationarity tests
   - Fills gaps
   > **Note**: Additional data cleaning, such as screening sensor flags, is recommended before using the partitioning class. While the package includes several quality control and assurance features, users are encouraged to perform their own tests as well. Refer to [Vickers and Mahrt, 1997](https://journals.ametsoc.org/view/journals/atot/14/3/1520-0426_1997_014_0512_qcafsp_2_0_co_2.xml) and [Zahn et al., 2016](http://article.sapub.org/10.5923.s.ajee.201601.20.html) for examples of additional tests. If you have suggestions on additional pre- and/or post-processing methods, please feel free to reach out or consider contributing to the code!

2. **Water-Use Efficiency Parameterizations**:
   Implements five parameterizations as described in [Zahn et al., 2021](https://www.sciencedirect.com/science/article/pii/S0168192321004767?via%3Dihub) "Direct Partitioning of Eddy-Covariance Water and Carbon Dioxide Fluxes into Ground and Plant Components".

3. **Partitioning Methods**:
   Implements five methods to partition eddy-covariance data and output T and E (W/m<sup>2</sup>) and P and R (mg_CO<sub>2</sub>/m<sup>2</sup>/s):
   
   1. **Conditional Eddy Covariance (CEC)**:
      [Zahn et al., 2021](https://www.sciencedirect.com/science/article/pii/S0168192321004767?via%3Dihub) "Direct Partitioning of Eddy-Covariance Water and Carbon Dioxide Fluxes into Ground and Plant Components".
   
   2. **Modified Relaxed Eddy Accumulation (MREA)**:
      [Thomas et al., 2008](https://www.sciencedirect.com/science/article/pii/S0168192308000737) "Estimating daytime subcanopy respiration from conditional sampling methods applied to multi-scalar high frequency turbulence time series".
   
   3. **Flux-Variance Similarity (FVS)**:
      [Scanlon and Sahu, 2008](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008WR006932) "Correlation-based flux partitioning of water vapor and carbon dioxide fluxes: Method simplification and estimation of canopy water use efficiency" and
      [Scanlon et al., 2019](https://www.sciencedirect.com/science/article/pii/S016819231930348X?via%3Dihub) "Correlation-based flux partitioning of water vapor and carbon dioxide fluxes: Method simplification and estimation of canopy water use efficiency".
   
   4. **Conditional Eddy Accumulation (CEA)**:
      [Zahn et al., 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JG008025) "Numerical Investigation of Observational Flux Partitioning Methods for Water Vapor and Carbon Dioxide".
   
   5. **Conditional Eddy Covariance with Water-Use Efficiency (CECw)**:
      [Zahn et al., 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JG008025) "Numerical Investigation of Observational Flux Partitioning Methods for Water Vapor and Carbon Dioxide".
      
4. **Method statistics**
    Implements statistics about the events sampled by each method on which the partitioning relies. [Thomas et al., 2008](https://www.sciencedirect.com/science/article/pii/S0168192308000737) "Estimating daytime subcanopy respiration from conditional sampling methods applied to multi-scalar high frequency turbulence time series"

5. **Easy data management**
    With one function call (`process`), large amount of data can be processed parallel with a variety of settings adjustable for each specific use case.




## Available Files

The following files are available in the repository:

- `src/partitioning/Partitioning.py`: Contains the `Partitioning` class with methods for quality control, data pre-processing, and five partitioning methods to separate ET and CO₂ fluxes into stomatal and non-stomatal components.
- `src/partitioning/auxfunctions.py`: Includes auxiliary functions for pre-processing and computing water-use efficiency, some of them adapted from [FluxPart](https://github.com/usda-ars-ussl/fluxpart).
- `src/wrapper.py`: Contains wrapper functions with a variety of settings to help loading, partition and save the data by one easy function call.

Old example files can be found in `tests/old_scripts_for_loading_and_processing`:
- `main.py`: Provides an example of how to use the script, with raw high-frequency data examples included in the `RawData30min` folder.
- `main_parallel.py`: Demonstrates how to run the script in parallel to process multiple files simultaneously.


---
## Setting up environment for developers

It is possible to install `PartitioningMethods` using pip in editable mode, where changes in the source code are directly mirrored in the installed package 
```bash
python -m pip install -e .
```

For developers it can be useful managing all the different environments with conda (this enables e.g. also using different python versions).
Using conda in the pre-devined environment I used:
 ```
 conda create -f conda/environment.yml
 conda activate PartMethods
 ```
In there you can develop and run from script, without installation.
For running from script without installation, the working directory just needs to be in the folder `src`.
Then its possible to just write  ```import partitioning``` in the python script.

If you want to install the package properly in there you first also need to install pip.
```
 conda install pip
 python -m pip install -e .
 ```

You should not have package version conflicts that way. If you do, use the environment version without any python versions specified and let conda do the rest:
```
 conda create -f conda/environment_noversions.yml
 conda activate PartMethods
 ```



---
## References for papers and datasets

- Zahn, E., Bou-Zeid, E., Good, S. P., Katul, G. G., Thomas, C. K., Ghannam, K., Smith, J. A., Chamecki, M., Dias, N. L., Fuentes, J. D., Alfieri, J. G., Kwon, H., Caylor, K. K., Gao, Z., Soderberg, K., Bambach, N. E., Hipps, L. E., Prueger, J. H., & Kustas, W. P. (2022). Direct partitioning of eddy-covariance water and carbon dioxide fluxes into ground and plant components. *Agricultural and Forest Meteorology, 315*, 108790. [https://doi.org/10.1016/j.agrformet.2021.108790](https://doi.org/10.1016/j.agrformet.2021.108790)

- Zahn, E., Ghannam, K., Chamecki, M., Moene, A. F., Kustas, W. P., Good, S. P., & Bou-Zeid, E. (2024). Numerical investigation of observational flux partitioning methods for water vapor and carbon dioxide. *Journal of Geophysical Research: Biogeosciences, 129*, e2024JG008025. [https://doi.org/10.1029/2024JG008025](https://doi.org/10.1029/2024JG008025)

- Zahn, E., & Bou-Zeid, E. (2024). Partitioning of water and CO2 fluxes at NEON sites into soil and plant components: A five-year dataset for spatial and temporal analysis. *Earth System Science Data Discussions* [preprint]. [https://doi.org/10.5194/essd-2024-272](https://doi.org/10.5194/essd-2024-272)

- Zahn, E., & Bou-Zeid, E. (2024). Partitioning of water and CO2 fluxes at NEON sites into soil and plant components: a five-year dataset for spatial and temporal analysis (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.12191876
