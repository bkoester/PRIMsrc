\name{PRIMsrc-package}

\alias{PRIMsrc}

\docType{package}

\title{Bump Hunting by Patient Rule Induction Method in Survival, Regression and Classification settings}

\description{
    Performs a unified treatment of Bump Hunting by Patient Rule Induction Method (PRIM) in Survival, Regression and Classification settings (SRC). 
    The method generates decision rules delineating a region in the predictor space, where the response is larger than its average over the entire space. 
    The region is shaped as a hyperdimensional box or hyperrectangle that is not necessarily contiguous. Assumptions are that the multivariate input 
    variables can be discrete or continuous and the univariate response variable can be discrete (Classification), continuous (Regression) 
    or a time-to event, possibly censored (Survival). It is intended to handle low and high-dimensional multivariate datasets, 
    including the situation where the number of covariates exceeds or dominates that of samples (\eqn{p > n} or \eqn{p \gg n} paradigm).
}

\details{
    The current version is a development release that only implements the case of a survival response. At this point, 
    survival bump hunting is also restricted to a directed peeling search of the first box covered by the recursive coverage (outer) loop of our 
    Patient Recursive Survival Peeling (PRSP) algorithm. New features will be added soon.
    
    The following describes the end-user functions that are needed to run a complete procedure.
    The other internal subroutines are not documented in the manual and are not to be called by the end-user at any time. 
    For computational efficiency, some end-user functions offer a parallelization option that is done by passing a few parameters 
    needed to configure a cluster. This is indicated by an asterisk (* = optionally involving cluster usage). 
    The R features are categorized as follows:

    \enumerate{

      \item{END-USER FUNCTION FOR PACKAGE NEWS} \cr
      \code{\link[PRIMsrc]{PRIMsrc.news}}
      \bold{Display the \pkg{PRIMsrc} Package News}\cr
            Function to display the log file \code{NEWS} of updates of the \pkg{PRIMsrc} package.\cr
            
      \item{END-USER S3-GENERIC FUNCTIONS FOR SUMMARY, DISPLAY, PLOT AND PREDICTION} \cr
      \code{\link[PRIMsrc]{summary}}
      \bold{Summary Function}\cr
            S3-generic summary function to summarize the main parameters used to generate the \code{PRSP} object.\cr
            
      \code{\link[PRIMsrc]{print}}
      \bold{Print Function}\cr
            S3-generic print function to display the cross-validated estimated values of the \code{PRSP} object.\cr

      \code{\link[PRIMsrc]{plot}}
      \bold{Two-Dimensional Visualization of Data Scatter and Box Vertices} \cr
            S3-generic plotting function for two-dimensional visualization of original or predicted data scatter 
            as well as cross-validated box vertices of a \code{PRSP} object. The scatter plot is for a given 
            peeling step of the peeling sequence and a given plane, both specified by the user.\cr

      \code{\link[PRIMsrc]{predict}}
      \bold{Predict Function}\cr
            S3-generic predict function to predict the box membership and box vertices 
            on an independent set.\cr
  
      \item{END-USER SURVIVAL BUMP HUNTING FUNCTIONS} \cr
      \code{\link[PRIMsrc]{cv.sbh}} (*)
      \bold{Cross-Validated Tuning of a Survival Bump Hunting Model} \cr
            First end-user function for tuning a cross-validated Survival Bump Hunting (SBH) model.
            Returns a cross-validated \code{CV} object, containing cross-validated estimates of optimal number of pre-selected variables 
            and model parameters. The function performs a single internal cross-validation procedure to simultaneously control model size 
            (#covariates) and model complexity (#peeling steps) before the model is fit. It does a univariate bump hunting variable 
            selection procedure, where model size and model complexity are simultaneously optimized using the cross-validation criterion 
            of choice: CER, LRT, or LHR (see companion paper below for details). This cross-validation procedure is carried 
            out separately from the main function \code{\link[PRIMsrc]{sbh}}. The returned S3-class \code{CV} object seves as input for 
            the main function \code{sbh} as well as the profiling plotting functions for the graphical display of profiling curves for model tuning
            (see \code{profile_plot} and \code{surface_plot} functions for more details). The function offers a number of options for the number 
            of cross-validation replicates to be perfomed: \eqn{B}; the type of cross-validation desired: \eqn{K}-fold (replicated)-averaged 
            or-combined, as well as the peeling and optimization critera chosen for model tuning and a few more parameters for the PRSP algorithm.
            The function takes advantage of the R package \pkg{parallel}, which allows users to create a cluster of workstations on a local 
            and/or remote machine(s), enabling scaling-up with the number of specified CPU cores and efficient parallel execution.\cr
      \code{\link[PRIMsrc]{sbh}} (*)
      \bold{Cross-Validated Fitting of a Survival Bump Hunting Model} \cr
            Second end-user function for fitting a cross-validated Survival Bump Hunting (SBH) model.
            Returns a cross-validated \code{PRSP} object, as generated by our Patient Recursive Survival Peeling or PRSP algorithm,
            containing cross-validated estimates of end-points statistics of interest. The function relies on a single internal 
            cross-validation procedure carried out by the main function \code{\link[PRIMsrc]{cv.sbh}} to simultaneously control 
            model size (#covariates) and model complexity (#peeling steps) before the model is fit. The returned S3-class \code{PRSP} 
            object contains cross-validated estimates of all the decision-rules of pre-selected covariates and all other statistical 
            quantities of interest at each iteration of the peeling sequence (inner loop of the PRSP algorithm). This enables the 
            graphical display of results of peeling trajectories, covariate traces and survival distributions 
            (see plotting functions for more details). The function takes advantage of the R package \pkg{parallel}, which allows
            users to create a cluster of workstations on a local and/or remote machine(s), enabling scaling-up with the number of 
            specified CPU cores and efficient parallel execution.\cr

      \item{END-USER PLOTTING FUNCTIONS FOR MODEL VALIDATION AND VISUALIZATION OF RESULTS} \cr
      \code{\link[PRIMsrc]{plot_profile}}
      \bold{One-Dimensional Visualization for Model Selection/Validation} \cr
            Function for plotting the cross-validated tuning profiles of a \code{CV} object.
            It uses the user's choice of statistics among the Log Hazard Ratio (LHR), Log-Rank Test (LRT) or Concordance Error Rate (CER) 
            as a function of the model tuning parameter, that is, the optimal number of peeling steps of the peeling sequence 
            (inner loop of our PRSP algorithm).\cr
      \code{\link[PRIMsrc]{plot_surface}}
      \bold{Two-Dimensional Visualization for Model Selection/Validation} \cr
            Function for plotting a perspective plot of the cross-validated tuning surface of a \code{CV} object. 
            It uses the user's choice of statistics among the Log Hazard Ratio (LHR), Log-Rank Test (LRT) or Concordance Error Rate (CER)
            as a function of the optimal choice/number of pre-selected variables and of the model tuning parameter, that is, and of number 
            of peeling steps of the peeling sequence (inner loop of our PRSP algorithm).\cr
      \code{\link[PRIMsrc]{plot_boxtraj}}
      \bold{Visualization of Peeling Trajectories/Profiles} \cr
            Function for plotting the cross-validated peeling trajectories/profiles of a \code{PRSP} object.
            Applies to the user-specified covariates among the pre-selected ones and all other statistical quantities of interest
            at each iteration of the peeling sequence (inner loop of our PRSP algorithm). \cr
      \code{\link[PRIMsrc]{plot_boxtrace}}
      \bold{Visualization of Covariates Traces} \cr
            Function for plotting the cross-validated covariates traces of a \code{PRSP} object.
            Plot the cross-validated modal trace curves of covariate importance and covariate usage of the user-specified 
            covariates among the pre-selected ones at each iteration of the peeling sequence (inner loop of our PRSP algorithm). \cr
      \code{\link[PRIMsrc]{plot_boxkm}}
      \bold{Visualization of Survival Distributions} \cr
            Function for plotting the cross-validated survival distributions of a \code{PRSP} object. 
            Plot the cross-validated Kaplan-Meir estimates of survival distributions for the highest risk (inbox) versus
            lower-risk (outbox) groups of samples at each iteration of the peeling sequence (inner loop of our PRSP algorithm). \cr

      \item{END-USER DATASETS} \cr
      \code{\link[PRIMsrc]{Synthetic.1}},
      \code{\link[PRIMsrc]{Synthetic.1b}}, 
      \code{\link[PRIMsrc]{Synthetic.2}},
      \code{\link[PRIMsrc]{Synthetic.3}}, 
      \code{\link[PRIMsrc]{Synthetic.4}}
      \bold{Five Datasets From Simulated Regression Survival Models} \cr
            Five datasets from simulated regression survival models #1-4 as described in Dazard et al. (2015) representing low- and high-dimensional situations, 
            and where regression parameters represent various types of relationship between survival times and covariates including saturated and noisy situations. 
            In three datasets where non-informative noisy covariates were used, these covariates were not part of the design matrix (models #2-3 and #4). 
            In one dataset, the signal is limited to a box-shaped region \eqn{R} of the predictor space (model #1b). 
            In the last dataset, the signal is limited to 10\% of the predictors in a \eqn{p > n} situation (model #4).
            Survival time was generated from an exponential model with with rate parameter \eqn{\lambda} (and mean \eqn{\frac{1}{\lambda}})
            according to a Cox-PH model with hazard exp(eta), where eta(.) is the regression function.
            Censoring indicator were generated from a uniform distribution on [0,3] (models #1-3) or [0,2] (model #4).
            In these synthetic datasets, all covariates are continuous, i.i.d. from a multivariate uniform distribution on [0,1] (models #1-3)
            or from a multivariate standard normal distribution (model #4).\cr
            
      \code{\link[PRIMsrc]{Real.1}}
      \bold{Clinical Dataset} \cr
            Publicly available HIV clinical data from the Women's Interagency HIV cohort Study (WIHS).
            Inclusion criteria of the study were that women at enrolment were (i) alive, (ii) HIV-1 infected, and 
            (iii) free of clinical AIDS symptoms. Women were followed until the first of the following occurred: 
            (i) treatment initiation (HAART), (ii) AIDS diagnosis, (iii) death, or administrative censoring. 
            The studied outcomes were the competing risks "AIDS/Death (before HAART)" and "Treatment Initiation (HAART)".
            However, here, for simplification purposes, only the first of the two competing events (i.e. the time to AIDS/Death), 
            was used in this dataset example. Likewise, the entire study enrolled 1164 women, but only the complete cases were used 
            in this dataset example for simplification. Variables included history of Injection Drug Use ("IDU") at enrollment,
            African American ethnicity ("Race"), age ("Age"), and baseline CD4 count ("CD4"). The question in this dataset example 
            was whether it is possible to achieve a prognostication of patients for AIDS and HAART.\cr

      \code{\link[PRIMsrc]{Real.2}}
      \bold{Genomic Dataset} \cr
            Publicly available lung cancer genomic data from the Chemores Cohort Study. This was an integrated study of mRNA, miRNA 
            and clinical variables to characterize the molecular distinctions between squamous cell carcinoma (SCC) 
            and adenocarcinoma (AC) in Non Small Cell Lung Cancer (NSCLC). Tissue samples were analysed from a cohort of 123 patients 
            who underwent complete surgical resection at the Institut Mutualiste Montsouris (Paris, France) between 30 January 2002 and 26 June 2006.
            In this genomic dataset, only the expression levels of Agilent miRNA probes (\eqn{p=939}) were included from the \eqn{n=123} samples of the Chemores cohort. 
            It represents a situation where the number of covariates dominates the number of complete observations, or \eqn{p >> n} case. \cr
    }
    Known Bugs/Problems : None at this time.
}

\author{
    \itemize{
        \item "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}
        \item "Michael Choe, M.D." \email{mjc206@case.edu}
        \item "Michael LeBlanc, Ph.D." \email{mleblanc@fhcrc.org}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}

    Acknowledgments: This project was partially funded by the National Institutes of Health
    NIH - National Cancer Institute (R01-CA160593) to J-E. Dazard and J.S. Rao.
}

\references{
    \itemize{
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{Cross-validation and Peeling Strategies for Survival Bump Hunting using Recursive Peeling Methods.}"
              Statistical Analysis and Data Mining (in press).
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2014).
              "\emph{Cross-Validation of Survival Bump Hunting by Recursive Peeling Methods.}"
              In JSM Proceedings, Survival Methods for Risk Estimation/Prediction Section. Boston, MA, USA.
              American Statistical Association IMS - JSM, p. 3366-3380.
        \item Dazard J-E., Choe M., LeBlanc M. and Rao J.S. (2015).
              "\emph{R package PRIMsrc: Bump Hunting by Patient Rule Induction Method for Survival, Regression and Classification.}"
              In JSM Proceedings, Statistical Programmers and Analysts Section. Seattle, WA, USA.
              American Statistical Association IMS - JSM, (in press).
        \item Dazard J-E. and J.S. Rao (2010).
              "\emph{Local Sparse Bump Hunting.}"
              J. Comp Graph. Statistics, 19(4):900-92.
    }
}

\keyword{Exploratory Survival/Risk Analysis}
\keyword{Survival/Risk Estimation & Prediction}
\keyword{Non-Parametric Method}
\keyword{Cross-Validation}
\keyword{Bump Hunting}
\keyword{Rule-Induction Method}

\seealso{
    \itemize{
        \item \code{makeCluster} (R package \pkg{parallel})
        \item \code{plot.survfit} (R package \pkg{survival})
        \item \code{glmnet} (R package \pkg{glmnet})
    }
}
