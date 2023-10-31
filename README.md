# Hirsch-et-al.-2023-test-of-the-diet-hypothesis-for-brain-expansion
This repository contains r code related to Hirsch et al 2023 test of the diet hypothesis for brain expansion

Original movement data are hosted on [Movebank](https://www.movebank.org/) (Processed data: Movebank ID 1120749252; Unprocessed data: Movebank ID 468460067)
 
 ## Instructions
There are various packages with complex dependencies used accross the scripts used in this study. In order to install every package, Windows users must first ensure that the latest version of [rtools](https://cran.r-project.org/bin/windows/Rtools/) is properly installed. Mac users must ensure that the latest version of [Xcode](https://developer.apple.com/xcode/) is properly installed.  This step cannot be skipped .

The following script installs the most recent versions of the required packages and their dependencies:

```
install.packages(c("devtools","miniCRAN","pacman"), dependencies = TRUE, type="source") 
devtools::install_github("wrathematics/getPass")
devtools::install_github("ctmm-initiative/ctmm")


#Check if required packages and their dependencies need installation or updates
list_of_required_packages <- c("plyr", "dplyr", "lubridate", "tidyr", "tidyverse", "data.table", "ggplot2", "ggpubr", "sp", 
"rgdal", "raster", "move", "lme4","emmeans", "StanHeaders", "rstan", "brms")

check_if_needs_install=as.character(miniCRAN::pkgDep(list_of_required_packages, suggests = TRUE, enhances = TRUE))
check_if_needs_update=as.character(pacman::p_update(FALSE))

new_packages <- check_if_needs_install[!(check_if_needs_install %in% installed.packages()[,"Package"])]
packages_to_update=check_if_needs_install[check_if_needs_install %in% check_if_needs_update]

packages_to_install=c(new_packages,packages_to_update)


install.packages(packages_to_install, type="source")
```

You will also need to install [cmdstan](https://mc-stan.org/docs/2_25/cmdstan-guide/cmdstan-installation.html) and [cmdstanr](https://mc-stan.org/cmdstanr/)


          
Be sure to change the working directory to an appropriate one on your system. Also be sure to modify where write.csv() saves all outputs. 

## Script metadata

```Reconstruct movement tracks to continuous time.R```: This script contains our implementation of path reconstruction to interpolate our data to continuous time
```FFT CTMM Distances.R```: This script implements the estimation of daily travel distance as described in Noonan et al. 2019
```Daily Path Length estimates.R```: This script estimates daily travel distance using the interpolated data as well as the original 4 min data and merges these estimates with the estimates from  ```FFT CTMM Distances.R```. The output of this script is used in ```Ben regression models.R```
```intertree stuff 2023.R``` This script is our implementation of estimating tree visitation and revisitation and uses the outputs of  ```Reconstruct movement tracks to continuous time.R```. 
```intertree stuff 2023_2_hour.R``` and ```intertree stuff 2023_4_hour.R``` are the same implementation as ```intertree stuff 2023.R``` but used on the first 2 and 4 hours of data respectively
```daily stuff for ben 2023 full day.R``` This script takes the output of ```intertree stuff 2023.R``` and calculates daily summaries of tree visitation for each animal.
```daily stuff for ben 2023 2 hour.R``` and ```daily stuff for ben 2023 4 hour.R``` are the same implementation as ```daily stuff for ben 2023 full day.R``` but used on the first 2 and 4 hours of data respectively
```Ben regression models.R``` takes the outputs of ```intertree stuff 2023.R```, ```intertree stuff 2023_2_hour.R```, ```intertree stuff 2023_4_hour.R```, ```daily stuff for ben 2023 full day.R```, ```daily stuff for ben 2023 2 hour.R``` and ```daily stuff for ben 2023 4 hour.R``` calculates efficiency metrics implements regression analyses.
```fft random movement comparison.R``` generates simulations of random movement models (brownian motion and OUF motion) and compares tree visitation rates of random movement to observed encounter rates. 


