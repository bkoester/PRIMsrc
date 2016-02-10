### General Remarks
<<<<<<< HEAD
=======

Grand-total number of CRAN downloads since initial release to CRAN (2015-07-28), 
as logged by [RStudio CRAN mirror](http://cran-logs.rstudio.com/):

![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc) 
 
Number of CRAN downloads in the last 30 days:

![](http://cranlogs.r-pkg.org/badges/PRIMsrc) 

Travis CI build result:

![](https://travis-ci.org/jedazard/PRIMsrc.svg)

===============
### Description

Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). 
The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. 
The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. 
>>>>>>> origin/devel

Grand-total number of CRAN downloads since initial release to CRAN (2015-07-28), 
as logged by [RStudio CRAN mirror](http://cran-logs.rstudio.com/):

[![](http://cranlogs.r-pkg.org/badges/grand-total/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

Number of CRAN downloads in the last month:

[![](http://cranlogs.r-pkg.org/badges/last-month/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)

Number of CRAN downloads in the last week:

[![](http://cranlogs.r-pkg.org/badges/last-week/PRIMsrc)](http://cran.rstudio.com/web/packages/PRIMsrc/index.html)


======================================
### Branch #1  (devel) - version 0.7.0

<<<<<<< HEAD
=======
======================================
### Branch #1  (devel) - version 0.7.0

>>>>>>> origin/devel
The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. 
Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. 
Specifically, it does a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion of choice: 
Concordance Error Rate (CER), Log-Rank Test (LRT), or Log-Hazard Ratio (LHR) (see companion paper below for details).

In addition, this cross-validation procedure is carried out in a cross-validation function called `cv.sbh()`, separately from the main function `sbh()`. 
Altogether, this allows a more rigorous treatment of model validation, a better control on the user-end and an improvement of the maintenance on the back-end. 
In the process, two S3-class objects are created instead of one: an additional S3-class object 'CV' is output by the cross-validation function `cv.sbh()` and used as input in the main function `sbh()`. 

<<<<<<< HEAD

=======
>>>>>>> origin/devel
===========
### License

PRIMsrc is Open Source / Free Software, available under the GNU General Public License, version 3. 
See details [here](https://github.com/jedazard/PRIMsrc/blob/devel/LICENSE).

<<<<<<< HEAD
=======
==============
### References

Open access to companion papers (accepted for publication):
>>>>>>> origin/devel

================
### Installation

* To install the most up-to-date version (0.7.0) of PRIMsrc from GitHub, using devtools:

`install.packages("devtools")`

`library("devtools")`

`devtools::install_github("jedazard/PRIMsrc")`


================
### Requirements

PRIMsrc 0.7.0 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-11-04 r69597) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:

[![Build Status](https://travis-ci.org/jedazard/PRIMsrc.png?branch=devel)](https://travis-ci.org/jedazard/PRIMsrc)


=========
### Usage

* To load the PRIMsrc library in an R session and start using it:

`library("PRIMsrc")`

* Check the package news with the R command:

`PRIMsrc.news()`

* Check on how to cite the package with the R command:

`citation("PRIMsrc")`

etc...


==============
### Wiki

<<<<<<< HEAD
See [Wiki](https://github.com/jedazard/PRIMsrc/wiki) page for Roadmap, Publications, Case Studies, Documentation and Manual, Examples and Support.
=======
- [JSM (2015)](https://www.amstat.org/membersonly/proceedings/2015/data/assets/pdf/233927.pdf). 
The ASA Proceedings of the annual Joint Statistical Meetings (Seattle, WA, USA).
>>>>>>> origin/devel

