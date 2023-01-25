[![DOI](https://zenodo.org/badge/593257800.svg)](https://zenodo.org/badge/latestdoi/593257800)
![publication](https://img.shields.io/badge/status-in%20review-orange)

# The influence of weather on the population dynamics of common mosquito vector species in the Canadian Prairies

*Submitted to Parasites & Vectors 2022*

Cole Baril<sup>1</sup>, Ben G. Pilling<sup>1</sup>, Milah J. Mikkelsen<sup>1</sup>, Jessica M. Sparrow<sup>1</sup>, Carlyn A. M Duncan<sup>1</sup>, Cody W. Koloski<sup>1</sup>, [Stefanie E. LaZerte](https://steffilazerte.ca)<sup>1,2</sup>, [Bryan J. Cassone](https://www.cassonelab.com/)<sup>1</sup>*

<sup>*</sup> Corresponding Author  
<sup>1</sup> Department of Biology, Brandon University, Brandon, MB R7A 6A9, Canada  
<sup>2</sup> Steffi LaZerte R Programming and Biological Consulting, Brandon, Manitoba  


> Please contact [Bryan Cassone](mailto:cassoneb@brandonu.ca) for questions regarding this paper and 
> Steffi LaZerte ([@steffilazerte](https://github.com/steffilazerte)) regarding 
> problems running these scripts.

# Contents

This repository contains the code used to load, prepare, explore and analyze 
mosquito counts for this manuscript.

## Format

The scripts are stored int the `Scripts` folder and can be run in sequence. 
Note that scripts are designed to produce formatted
RMarkdown reports, so can be run interactively (by hand) or with the `rmarkdown::render()` 
code contained at the top of each script. If running this code, an html report
of the script will be produced in the `Results` folder. 
Data produced will be saved to a `Data` folder. 
The original data is stored in the `Data/Raw` folder.
The final dataset, `mosquitoes.csv` is stored in the `Data/Datasets` folder.

To improve reproducibility, [`renv`](https://rstudio.github.io/renv) is used to 
set the package versions used. To restore this set of packages, open this project
as an Rstudio project, and use the `renv::restore()` function (listed but commented
out at the top of each script).

