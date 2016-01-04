### General Remarks

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


Assumptions are that the multivariate input covariates can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event, possibly censored (Survival).
It is intended to handle low and high-dimensional multivariate datasets, including the situation where the number of covariates exceeds or dominates that of samples (p > n or p >> n paradigm).

========================================
### Branch #2  (unified) - version 1.0.0

This second branch (unified) will host the future complete version of the code (version 1.0.0), including undirected peeling search by Patient Rule Induction Method (PRIM) 
that will allow the unified treatment of bump hunting for every type of common response: Survival, Regression and Classification (SRC).

===========
### License

PRIMsrc is Open Source / Free Software, available under the GNU General Public License, version 3. 
See details [here](https://github.com/jedazard/PRIMsrc/blob/unified/LICENSE).

==============
### References

Open access to companion papers (accepted for publication):

- [Statistical Analysis and Data Mining (2015-12-09)](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1932-1872). 
The American Statistical Association (ASA) Affiliated Data Science Journal (to appear).

- [arXiv v1:v8 (2015-01-16 : 2015-11-20)](http://arxiv.org/abs/1501.03856). 
The Cornell University Library Archives.

- [JSM (2015)](https://www.amstat.org/membersonly/proceedings/2015/data/assets/pdf/233927.pdf). 
The ASA Proceedings of the annual Joint Statistical Meetings (Seattle, WA, USA).

- [JSM (2014)](https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf). 
The ASA Proceedings of the annual Joint Statistical Meetings (Boston, MA, USA).
