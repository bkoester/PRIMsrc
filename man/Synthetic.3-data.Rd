\name{Synthetic.3}

\docType{data}

\alias{Synthetic.3}

\title{Synthetic Dataset #3: \eqn{p < n} case}

\description{
  Dataset from simulated regression survival model #3 as described in Dazard et al. (2015). 
  Here, the regression function does not include any of the predictors. 
  This means that none of the covariates is informative (noisy), and are not part of the design matrix. 
  Survival time was generated from an exponential model with rate parameter \eqn{\lambda} (and mean \eqn{\frac{1}{\lambda}})
  according to a Cox-PH model with hazard exp(eta), where eta(.) is the regression function.
  Censoring indicator were generated from a uniform distribution on [0, 3].
  In this synthetic example, all covariates are continuous, i.i.d. from a multivariate uniform distribution on [0, 1].
}

\usage{Synthetic.3}

\format{
  Each dataset consists of a \code{numeric} \code{matrix} containing \eqn{n=250} observations (samples) 
  by rows and \eqn{p=3} variables by columns, not including the censoring indicator and (censored) time-to-event variables.
  It comes as a compressed Rda data file.
}

\source{See simulated survival model #3 in Dazard et al., 2015.}

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

\keyword{datasets}

