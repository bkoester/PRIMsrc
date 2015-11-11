============
Description:
============
Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. Assumptions are that the multivariate input covariates can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) or a time-to event, possibly censored (Survival). It is intended to handle low and high-dimensional multivariate datasets, including the situation where the number of covariates exceeds or dominates that of samples (p > n or p >> n paradigm).

===================================
Branch #1  (devel) - version 0.7.0:
===================================
The first branch (devel) hosts a development version of the code (version 0.7.0) that is more rigorous and modular. Here, a single internal cross-validation procedure is carried out to simultaneously control model size (#covariates) and model complexity (#peeling steps) before the model is fit. Specifically, it includes a univariate bump hunting variable selection procedure, where model size and model complexity are simultaneously optimized by cross-validation of the cross-validation criterion of choice: CER, LRT, or LHR (see companion paper below for details).

In addition, this cross-validation procedure is carried out separately of the main function 'sbh()' in a cross-validation function called 'cv.sbh()'. Altogether, this allows a more rigorous treatment of model validation, a better control on the user-end and an improvement of the maintenance on the back-end. In the process, two S3-class objects are created instead of one: an additional S3-class object 'CV' is output by the cross-validation function cv.sbh() and used as input in the main function 'sbh()'. 

========
License:
========
PRIMsrc is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

===========
References:
===========
The companion papers (accepted and submitted for publication) can be accessed here:

- ASA-IMS JSM Proceedings (2014): 

https://www.amstat.org/membersonly/proceedings/2014/data/assets/pdf/312982_90342.pdf

- Archives arXiv:

http://arxiv.org/abs/1501.03856.

- Statistical Analysis and Data Mining. The ASA Data Science Journal (to appear):

http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1932-1872
