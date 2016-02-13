### General Remarks

CRAN downloads since initial release to CRAN (2015-07-28):
[![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)
as tracked by [RStudio CRAN mirror](http://cran-logs.rstudio.com/)

CRAN downloads in the last month:
[![](http://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

CRAN downloads in the last week:
[![](http://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)


======================================
### Branch #1  (devel) - version 0.7.0

The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. 
Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. 
Specifically, it does a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion of choice: 
Concordance Error Rate (CER), Log-Rank Test (LRT), or Log-Hazard Ratio (LHR) (see companion paper below for details).

In addition, this cross-validation procedure is carried out in a cross-validation function called `cv.sbh()`, separately from the main function `sbh()`. 
Altogether, this allows a more rigorous treatment of model validation, a better control on the user-end and an improvement of the maintenance on the back-end. 
In the process, two S3-class objects are created instead of one: an additional S3-class object 'CV' is output by the cross-validation function `cv.sbh()` and used as input in the main function `sbh()`. 


===========
### License

PRIMsrc is open source / free software, licensed under the GNU General Public License, version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


================
### Installation

* To install the most up-to-date development version (0.7.0) of PRIMsrc from GitHub, using devtools:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jedazard/PRIMsrc")
```


================
### Requirements

PRIMsrc 0.7.0 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-11-04 r69597) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:

[![Build Status](https://travis-ci.org/jedazard/PRIMsrc.png?branch=devel)](https://travis-ci.org/jedazard/PRIMsrc)


=========
### Usage

* To load the PRIMsrc library in an R session and start using it:

```{r}
library("PRIMsrc")
```

* Check the package news with the R command:

```{r}
PRIMsrc.news()
```

* Check on how to cite the package with the R command:

```{r}
citation("PRIMsrc")
```

etc...


==========================
### Website - Wiki

- See [Website](http://jedazard.github.io/PRIMsrc/) 
- See [Wiki](https://github.com/jedazard/PRIMsrc/wiki) for Roadmap, Publications, Case Studies, Documentation and Manual, Examples and Support.
