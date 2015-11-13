\name{cv.sbh}

\alias{cv.sbh}

\title{Cross-Validated Tuning of a Survival Bump Hunting Model}

\description{
   First end-user function for tuning a cross-validated Survival Bump Hunting (SBH) model.
   Returns a cross-validated \code{CV} object, containing cross-validated estimates 
   of optimal number of pre-selected variables and model parameters.
}

\usage{
  cv.sbh(dataset, 
         B = 10, K = 5, 
         cvtype = c("combined", "averaged", "none", NULL), 
         cvcriterion = c("lrt", "cer", "lhr", NULL),
         conservative = c("most", "medium", "least"),
         arg = "beta=0.05,alpha=0.05,minn=5,peelcriterion=\"lr\"",
         fdr = NULL, thr = NULL, 
         parallel = FALSE, conf = NULL, seed = NULL)
}

\arguments{
  \item{dataset}{
    \code{data.frame} or \code{numeric} \code{matrix} of input dataset containing the observed survival and status indicator variables 
    in the first two columns, respectively, and all the covariates thereafter. If a \code{data.frame} is provided, it will be coerced 
    to a \code{numeric} \code{matrix}. Discrete (or nominal) covariates should be made (or re-arranged into) ordinal variables.
  }
  \item{B}{
    Positive \code{integer} scalar of the number of replications of the cross-validation procedure. Defaults to 10.
  }
  \item{K}{
    \code{Integer} giving the number of folds (partitions) into which the observations should be randomly split 
    for the cross-validation procedure. Setting K also specifies the type of cross-validation to be done:
    \itemize{
        \item{\code{K}} = 1 carries no cross-validation out.
        \item{\code{K}} in \{2,...,\eqn{n-1}\} carries out eqn{K}-fold cross-validation.
        \item{\code{K}} = \eqn{n} carries out leave-one-out cross-validation.
    }
  } 
  \item{cvtype}{
    \code{Character} \code{vector} describing the cross-validation technique in \{"combined", "averaged", "none", NULL\}.
    If \code{NULL}, automatically reset to "none".
  }
  \item{cvcriterion}{
    \code{character} \code{vector} describing the optimization criterion in \{"lrt", "lhr", "cer", NULL\}.
    If \code{NULL}, automatically reset to "none".
  }
  \item{conservative}{
    \code{character} \code{vector} describing the degree of conservativeness in \{"most", "medium", "least"\}
    to be used in variable pre-selection.
  }
  \item{arg}{
    \code{Character} \code{vector} describing the PRSP parameters:
    \itemize{
        \item{\code{alpha}} = fraction to peel off at each step. Defaults to 0.05.
        \item{\code{beta}} = minimum support size resulting from the peeling sequence. Defaults to 0.05.
        \item{\code{minn}} = minimum number of observation that we want to be able to detect in a box. Defaults to 5.
        \item{\code{peelcriterion}} in \{"hr" for Log-Hazard Ratio (LHR), 
                                         "lr" for Log-Rank Test (LRT), 
                                         "ch" for Cumulative Hazard Summary (CHS)\}. 
                    Defaults to "lr".
    }
    Note that the parameters in \code{arg} come as a string of charaters between double quotes, 
    where all parameter evaluations are separated by comas (see example).
  } 
  \item{fdr}{
    \code{Numeric} scalar of the FDR level at which we want to pre-select the variables. Defaults to NULL.
  }
  \item{thr}{
    \code{Numeric} scalar of the threshold number of pre-selected the variables. Defaults to NULL.
  }
  \item{parallel}{
    \code{Logical}. Is parallel computing to be performed? Optional. Defaults to \code{FALSE}.
  }
  \item{conf}{
    \code{List} of parameters for cluster configuration.
    Inputs for R package \pkg{parallel} function \code{makeCluster} (R package \pkg{parallel}) for cluster setup.
    Optional, defaults to \code{NULL}. See details for usage.
  }
  \item{seed}{
    Positive \code{integer} scalar of the user seed to reproduce the results.
 }
}

\details{
  \code{cv.sbh} performs a single internal cross-validation procedure to simultaneously control model size (#covariates) 
  and model complexity (#peeling steps) before the model is fit. It does a univariate bump hunting variable selection 
  procedure, where model size and model complexity are simultaneously optimized using the cross-validation 
  criterion of choice: CER, LRT, or LHR (see companion paper below for details). This cross-validation procedure 
  is carried out separately from the main function \code{\link[PRIMsrc]{sbh}}. The returned S3-class \code{CV} object 
  seves as input for the main function \code{sbh} as well as the profiling plotting functions for the graphical display 
  of profiling curves for model tuning (see \code{profile_plot} and \code{surface_plot} functions for more details).

  The function offers a number of options for the number of cross-validation replicates to be perfomed: \eqn{B}; 
  the type of cross-validation desired: \eqn{K}-fold (replicated)-averaged or-combined, as well as the peeling and optimization 
  critera chosen for model tuning and a few more parameters for the PRSP algorithm.
    
  The function takes advantage of the R package \pkg{parallel}, which allows users to create a cluster of workstations on a local
  and/or remote machine(s), enabling scaling-up with the number of specified CPU cores and efficient parallel execution.
  
  In the case of large (\eqn{p > n}) or very large (\eqn{p >> n}) datasets, it is required to use the parallelization option.
  To run a parallel session (and parallel RNG) of the PRIMsrc procedures (\code{parallel}=\code{TRUE}), argument \code{conf}
  is to be specified (i.e. non \code{NULL}). It must list the specifications of the folowing parameters for cluster configuration:
  "names", "cpus", "type", "homo", "verbose", "outfile". These match the arguments described in function \code{makeCluster}
  of the R package \pkg{parallel}. All fields are required to properly configure the cluster, except for "names" and "cpus",
  which are the values used alternatively in the case of a cluster of type "SOCK" (socket), or in the case of a cluster
  of type other than "SOCK" (socket), respectively. See examples below.
  \itemize{
        \item "names": \code{names} : \code{character} \code{vector} specifying the host names on which to run the job.
              Could default to a unique local machine, in which case, one may use the unique host name "localhost".
              Each host name can potentially be repeated to the number of CPU cores available on the corresponding machine.
        \item "cpus": \code{spec} : \code{integer} scalar specifying the total number of CPU cores to be used
              across the network of available nodes, counting the workernodes and masternode.
        \item "type": \code{type} : \code{character} \code{vector} specifying the cluster type ("SOCK", "PVM", "MPI").
        \item "homo": \code{homogeneous} : \code{logical} scalar to be set to \code{FALSE} for inhomogeneous clusters.
        \item "verbose": \code{verbose} : \code{logical} scalar to be set to \code{FALSE} for quiet mode.
        \item "outfile": \code{outfile} : \code{character} \code{vector} of the output log file name for the workernodes.
  }

  Note that argument \code{B} is internally reset to \code{conf$cpus}*\code{ceiling}(\code{B}/\code{conf$cpus}) in case the 
  parallelization is used (i.e. \code{conf} is non \code{NULL}), where \code{conf$cpus} denotes the total number of CPUs to be 
  used (see above).
  
  The actual creation of the cluster, its initialization, and closing are all done internally.
  In addition, when random number generation is needed, the creation of separate streams of parallel RNG per node
  is done internally by distributing the stream states to the nodes (For more details see function \code{makeCluster}
  (R package \pkg{parallel}) and/or \url{http://www.stat.uiowa.edu/~luke/R/cluster/cluster.html}.
  
  The use of a seed allows to reproduce the results within the same type of session: the same seed will reproduce the same results 
  within a non-parallel session or within a parallel session, but it will not necessarily give the exact same results (up to sampling variability)
  between a non-parallelized and parallelized session due to the difference of management of the seed between the two (see parallel RNG and 
  value of retuned seed below).
}

\value{
  Object of \code{class} \code{PRSP} (Patient Recursive Survival Peeling)
  \code{List} containing the following 19 fields:
    \item{x}{\code{numeric} \code{matrix} of original dataset.}
    \item{times}{\code{numeric} \code{vector} of observed failure / survival times.}
    \item{status}{\code{numeric} \code{vector} of observed event indicator in \{1,0\}.}
    \item{B}{positive \code{integer} of the number of replications used in the cross-validation procedure.}
    \item{K}{positive \code{integer} of the number of folds used in the cross-validation procedure.}
    \item{cvtype}{\code{character} \code{vector} of the cross-validation technique used.}
    \item{cvcriterion}{\code{character} \code{vector} of optimization criterion used.}
    \item{conservative}{\code{character} \code{vector} degree of conservativeness used.}
    \item{arg}{\code{character} \code{vector} of the parameters used.}
    \item{fdr}{\code{character} \code{vector} of the \code{fdr} parameter used.}
    \item{thr}{\code{character} \code{vector} of the \code{thr} parameter used.}
    \item{cvfit}{\code{List} with 7 fields of cross-validated estimates:
        \itemize{
            \item{cv.maxsteps}: \code{numeric} scalar of maximal number of peeling steps over the replicates.
            \item{cv.nsteps}: \code{numeric} scalar of optimal number of peeling steps according to the optimization criterion.
            \item{cv.model}: \code{numeric} scalar of coordinate of optimal model.
            \item{cv.nmodel}: \code{numeric} scalar of number of retained models or numbers of pre-selected variables.
            \item{cv.sel}{\code{numeric} \code{vector} of pre-selected covariates in reference to original index.}
            \item{cv.sign}{\code{numeric} \code{vector} in \{-1,+1\} of directions of peeling for all pre-selected covariates.} 
            \item{cv.profile}{\code{List} of 3 \code{numeric} objects: 
                  - an \code{array} of matrices, one for each cross-validation replicate (in the third dimension), 
                    of the cross-validated statistics used in the optimization criterion (set by user) as a function of 
                    the number of pre-selected variales and the number of peeling steps.  
                  - a \code{matrix} of pointwise mean values of the above array over the cross-validation replicates.
                  - a \code{matrix} of pointwise standard error values of the above array over the cross-validation replicates.}
        }
    \item{success}{\code{logical} scalar of the returned flag for pursuing or not a SBH model fitting by the \code{sbh} function.}
    \item{seed}{User seed(s) used:
        \code{integer} of a single value, if parallelization is used
        \code{integer} \code{vector} of values, one for each replication, if parallelization is not used.}
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

\note{
    Unique end-user function for fitting the Survival Bump Hunting model.
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
        \item \code{cv.glmnet} (R package \pkg{glmnet})
        \item \code{glmnet} (R package \pkg{glmnet})
    }
}

\examples{
#===================================================
# Loading the library and its dependencies
#===================================================
library("PRIMsrc")

#===================================================
# Package news
# Package citation
#===================================================
PRIMsrc.news()
citation("PRIMsrc")
    
#===================================================
# Demo with a synthetic dataset
# Use help for descriptions
#===================================================
data("Synthetic.1", package="PRIMsrc")
?Synthetic.1

#===================================================
# Simulated dataset #1 (n=250, p=3)
# Replicated Combined Cross-Validation (RCCV)
# Peeling criterion = LRT
# Optimization criterion = LRT
# Without parallelization
# Without computation of permutation p-values
#===================================================
CVCOMBREP.CV <- cv.sbh(dataset = Synthetic.1, 
                       B = 1, K = 5, 
                       cvtype = "combined", 
                       cvcriterion = "lrt",
                       conservative = "least",
                       arg = "beta=0.05,
                              alpha=0.05,
                              minn=5,
                              peelcriterion=\"lr\"",
                       fdr = NULL, thr = NULL,
                       parallel = FALSE, conf = NULL, seed = 123)

\dontrun{
    #===================================================
    # Examples of parallel backend parametrization 
    #===================================================
    # Example #1 - 1-Quad (4-core double threaded) PC 
    # Running WINDOWS
    # With SOCKET communication
    #===================================================
    if (.Platform$OS.type == "windows") {
        cpus <- detectCores()
        conf <- list("names" = rep("localhost", cpus),
                     "cpus" = cpus,
                     "type" = "SOCK",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    #===================================================
    # Example #2 - 1 master node + 3 worker nodes cluster
    # All nodes equipped with identical setups and multicores
    # Running LINUX
    # With SOCKET communication
    #===================================================
    if (.Platform$OS.type == "unix") {
        masterhost <- Sys.getenv("HOSTNAME")
        slavehosts <- c("compute-0-0", "compute-0-1", "compute-0-2")
        nodes <- length(slavehosts) + 1
        cpus <- 8
        conf <- list("names" = c(rep(masterhost, cpus),
                                 rep(slavehosts, cpus)),
                     "cpus" = nodes * cpus,
                     "type" = "SOCK",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    #===================================================
    # Example #3 - Multinode multicore per node cluster
    # Running LINUX 
    # with MPI communication
    # Here, a file named ".nodes" (e.g. in the home directory)
    # contains the list of nodes of the cluster
    #===================================================
    if (.Platform$OS.type == "unix") {
        hosts <- scan(file=paste(Sys.getenv("HOME"), "/.nodes", sep=""), 
                      what="", 
                      sep="\n")
        hostnames <- unique(hosts)
        nodes <- length(hostnames)
        cpus <-  length(hosts)/length(hostnames)
        conf <- list("cpus" = nodes * cpus,
                     "type" = "MPI",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    
    #===================================================
    # Simulated dataset #1 (n=250, p=3)
    # Replicated Combined Cross-Validation (RCCV)
    # Peeling criterion = LRT
    # Optimization criterion = LRT
    # With parallelization
    # With computation of permutation p-values
    #===================================================
    CVCOMBREP.CV <- cv.sbh(dataset = Synthetic.1, 
                       B = 1, K = 5, 
                       cvtype = "combined", 
                       cvcriterion = "lrt",
                       conservative = "least",
                       arg = "beta=0.05,
                              alpha=0.05,
                              minn=5,
                              peelcriterion=\"lr\"",
                       fdr = NULL, thr = NULL,
                       parallel = TRUE, conf = conf, seed = 123)

}
}
