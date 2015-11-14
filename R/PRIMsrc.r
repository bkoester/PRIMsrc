##########################################################################################################################################
# PRIMsrc
##########################################################################################################################################

##########################################################################################################################################
# 1. END-USER SURVIVAL BUMP HUNTING FUNCTIONS
##########################################################################################################################################

##########################################################################################################################################
################
# Usage         :
################
#                   cv.sbh(dataset,
#                          B=10, K=5,
#                          cvtype=c("combined", "averaged", "none", NULL),
#                          cvcriterion=c("lrt", "cer", "lhr", NULL),
#                          conservative=c("most", "medium", "least"),
#                          arg="beta=0.05,alpha=0.05,minn=5,peelcriterion=\"lr\"",
#                          fdr=NULL,
#                          thr=NULL,
#                          parallel=FALSE, conf=NULL, seed=NULL)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.sbh <- function(dataset,
                   B=10, K=5,
                   cvtype=c("combined", "averaged", "none", NULL),
                   cvcriterion=c("lrt", "cer", "lhr", NULL),
                   conservative=c("most", "medium", "least"),
                   arg="beta=0.05,alpha=0.05,minn=5,peelcriterion=\"lr\"",
                   fdr=NULL,
                   thr=NULL,
                   parallel=FALSE, conf=NULL, seed=NULL) {

  # Parsing and evaluating parameters
  alpha <- NULL
  beta <- NULL
  minn <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  # Checking matching parameters
  cvtype <- match.arg(arg=cvtype, choices=c(NULL, "combined", "averaged", "none"), several.ok=FALSE)
  cvcriterion <- match.arg(arg=cvcriterion, choices=c(NULL, "lrt", "cer", "lhr"), several.ok=FALSE)
  conservative <- match.arg(arg=conservative, choices=c("most", "medium", "least"), several.ok=FALSE)

  # Checking inputs
  if (missing(dataset)) {
    stop("\nNo dataset provided !\n\n")
  } else {
    cat("\nSurvival dataset provided.\n\n")
    digits <- getOption("digits")
    if (!(is.data.frame(dataset))) {
      dataset <- as.data.frame(dataset)
    }
    x <- as.matrix(dataset[ ,-c(1,2), drop=FALSE])
    mode(x) <- "numeric"
    n <- nrow(x)
    p <- ncol(x)
    mini <- max(minn, ceiling(beta*n))
    if (n < mini) {
        stop("Based on the input parameters, the number of data points must be greater than the threshold of ", mini, " points.\n")
    }
    times <- dataset$stime
    status <- dataset$status
    times[times <= 0] <- 10^(-digits)
    if (is.null(colnames(x))) {
        colnames(x) <- paste("X", 1:p, sep="")
    }
  }

  if ((is.null(cvtype)) || (cvtype == "none") || (is.null(cvcriterion)) || (cvcriterion == "none")) {
    cvtype <- "none"
    cvcriterion <- "none"
    fdr <- NULL
    thr <- NULL
    B <- 1
    K <- 1
    conf <- NULL
    parallel <- FALSE
    cat("No cross-validation requested. No replication will be performed. No need of parallelization. \n")
  } else {
    if (B > 1) {
        if (parallel) {
            cat("Requested parallel replicated ", K, "-fold cross-validation with ", conf$cpus*ceiling(B/conf$cpus), " replications \n", sep="")
        } else {
            cat("Requested replicated ", K, "-fold cross-validation with ", B, " replications \n", sep="")
        }
    } else {
        cat("Requested single ", K, "-fold cross-validation without replications \n", sep="")
    }
  }
  cat("Variable pre-selection conservativeness: ", disp(x=conservative), "\n")
  cat("Cross-validation technique: ", disp(x=cvtype), "\n")
  cat("Cross-validation criterion: ", disp(x=cvcriterion), "\n")
  cat("Peeling criterion: ", disp(x=peelcriterion), "\n")
  cat("Parallelization:", parallel, "\n")
  cat("\n")

  if (is.null(fdr) && is.null(thr)) {
    M <- 1
  } else {
    if (is.null(thr)) {
        M <- length(fdr)
    } else if (is.null(fdr)) {
        M <- length(thr)
    }
  }

  # CV of Survival Bump Hunting model
  cat("Cross-Validation and variable selection of Survival Bump Hunting model using the PRSP algorithm ... \n")
  if (!parallel) {
      if (is.null(seed)) {
          seed <- runif(n=B, min=1, max=2) * 10^(digits-2)
      } else {
          seed <- (0:(B-1)) + seed
      }
      CV.peel.rep.obj <- cv.tune.rep(x=x, times=times, status=status,
                                     B=B, K=K, arg=arg,
                                     fdr=fdr,
                                     thr=thr,
                                     cvtype=cvtype,
                                     cvcriterion=cvcriterion,
                                     parallel=parallel, seed=seed)
  } else {
     if (conf$type == "SOCK") {
          cl <- makeCluster(spec=conf$names,
                            type=conf$type,
                            homogeneous=conf$homo,
                            outfile=conf$outfile,
                            verbose=conf$verbose)
      } else {
          cl <- makeCluster(spec=conf$cpus,
                            type=conf$type,
                            homogeneous=conf$homo,
                            outfile=conf$outfile,
                            verbose=conf$verbose)
      }
      clusterSetRNGStream(cl=cl, iseed=seed)
      a <- ceiling(B/conf$cpus)
      B <- a*conf$cpus
######################################################################################################
clusterEvalQ(cl=cl, expr=library("survival"))
clusterEvalQ(cl=cl, expr=library("Hmisc"))
clusterExport(cl=cl,
              varlist=list("cv.tune.rep",
                           "cv.ave.tune", "cv.comb.tune",
                           "cv.folds", "cv.tune", "prsp",
                           "updatecut", "endpoints",
                           "is.empty", "cbindlist",
                           "list2mat", "list2array",
                           "lapply.mat", "lapply.array"),
              envir=.GlobalEnv)
######################################################################################################
      obj.cl <- clusterCall(cl=cl, fun=cv.tune.rep,
                            x=x, times=times, status=status,
                            B=a, K=K, arg=arg,
                            fdr=fdr,
                            thr=thr,
                            cvtype=cvtype,
                            cvcriterion=cvcriterion,
                            parallel=parallel, seed=NULL)
      stopCluster(cl)
      CV.peel.rep.obj <- list("cv.maxsteps"=vector(mode="list", length=B),
                              "cv.sel"=vector(mode="list", length=B),
                              "cv.sign"=vector(mode="list", length=B),
                              "cv.profile"=vector(mode="list", length=B))
      for (b in 1:conf$cpus) {
          CV.peel.rep.obj$cv.maxsteps[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.maxsteps
          CV.peel.rep.obj$cv.sel[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.sel
          CV.peel.rep.obj$cv.sign[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.sign
          CV.peel.rep.obj$cv.profile[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.profile
      }
      CV.peel.rep.obj$success <- obj.cl[[1]]$success
  }

  # Collect the peeling statistics for each step from all the replicates
  CV.maxsteps.rep <- CV.peel.rep.obj$cv.maxsteps
  CV.sel.rep <- CV.peel.rep.obj$cv.sel
  CV.sign.rep <- CV.peel.rep.obj$cv.sign
  CV.profile.rep <- CV.peel.rep.obj$cv.profile
  success <- CV.peel.rep.obj$success

  if (success) {

        cat("Success! ", B, " (replicated) cross-validation(s) has(ve) completed \n", sep="")

        # Cross-validated maximum peeling length, thresholded by minimal box support, from all replicates
        CV.maxsteps <- ceiling(mean(unlist(CV.maxsteps.rep)))

        # CV profiles with mean and standard error from all replicates
        cat("Generating cross-validated profiles and optimal peeling length ...\n")
        if (cvtype == "none") {
            CV.profile.array <- NULL
            CV.mean.profile <- NULL
            CV.se.profile <- NULL
        } else if ((cvtype == "averaged") || (cvtype == "combined")) {
            CV.profile.array <- list2array(list=CV.profile.rep, fill=NA, coltrunc=CV.maxsteps)
            dimnames(CV.profile.array) <- list(paste("model", 1:M, sep=""), paste("step", 0:(CV.maxsteps-1), sep=""), 1:B)
            CV.mean.profile <- apply(X=CV.profile.array, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
            CV.se.profile <- apply(X=CV.profile.array, MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
            dimnames(CV.mean.profile) <- list(paste("model", 1:M, sep=""), paste("step", 0:(CV.maxsteps-1), sep=""))
            dimnames(CV.se.profile) <- list(paste("model", 1:M, sep=""), paste("step", 0:(CV.maxsteps-1), sep=""))
        } else {
            stop("Invalid CV type option \n")
        }
        CV.profile <- list("profiles"=CV.profile.array,
                           "mean"    =CV.mean.profile,
                           "se"      =CV.se.profile)

        # Cross-validated optimal peeling length over the replicates
        if (cvtype == "none") {
            CV.model <- 1
            CV.nsteps <- CV.maxsteps
        } else {
            z <- CV.profile$mean
            if ((cvcriterion == "lhr") || (cvcriterion == "lrt")) {
                w <- which(x=(z == as.numeric(z)[which.max(z)]), arr.ind=TRUE, useNames=TRUE)
            } else if (cvcriterion == "cer") {
                w <- which(x=(z == as.numeric(z)[which.min(z)]), arr.ind=TRUE, useNames=TRUE)
            }
            CV.model <- w[1,1]
            CV.nsteps <- w[1,2]
        }

        # Common selected covariates over the replicates
        cat("Generating common selected covariates over the replicates ...\n")
        CV.sel <- vector(mode="list", length=M)
        names(CV.sel) <- paste("model", 1:M, sep="")
        varsel <- unique(unlist(CV.sel.rep))
        maxp <- length(varsel)
        wbm <- vector(mode="list", length=B)
        CV.sel.B <- array(dim=c(M,maxp,B), data=NA, dimnames=list(paste("model", 1:M, sep=""), colnames(x)[varsel], 1:B))
        for (b in 1:B) {
            wbm[[b]] <- vector(mode="list", length=M)
            for (m in 1:M) {
                wbm[[b]][[m]] <- pmatch(x=CV.sel.rep[[b]][[m]], table=varsel, nomatch=NA, duplicates.ok=FALSE)
                CV.sel.B[m,wbm[[b]][[m]],b] <- CV.sel.rep[[b]][[m]]
            }
        }
        for (m in 1:M) {
            varfreq <- table(CV.sel.B[m,,], useNA="no")
            if (conservative == "most") {
                CV.sel[[m]] <- as.numeric(names(varfreq[varfreq == B]))              # most conservative (any covariate present in all the replications)
            } else if (conservative == "medium") {
                CV.sel[[m]] <- as.numeric(names(varfreq[varfreq >= ceiling(B/2)]))   # medium conservative (any covariate present in at least half of the replications)
            } else if (conservative == "least") {
                CV.sel[[m]] <- as.numeric(names(varfreq))                            # least conservative (any covariate present in at least one of the replications)
            }
            names(CV.sel[[m]]) <- colnames(x)[CV.sel[[m]]]
        }

        # Modal sign values of selected covariates over the replicates (by step)
        cat("Generating modal sign values of selected covariates over the replicates ...\n")
        CV.sign <- vector(mode="list", length=M)
        names(CV.sign) <- paste("model", 1:M, sep="")
        CV.sign.B <- array(dim=c(M,maxp,B), data=NA, dimnames=list(paste("model", 1:M, sep=""), colnames(x)[varsel], 1:B))
        for (b in 1:B) {
            for (m in 1:M) {
                CV.sign.B[m,wbm[[b]][[m]],b] <- CV.sign.rep[[b]][[m]]
            }
        }
        signmode <- apply(X=CV.sign.B,
                          MARGIN=c(1,2),
                          FUN=function(x){if (all(is.na(x)))
                                            return(NA)
                                          else
                                            return(as.numeric(names(which.max(table(x, useNA="no")))))
                                         }
                         )
        for (m in 1:M) {
            CV.sign[[m]] <- signmode[m,names(CV.sel[[m]])]
        }

  } else {

        cat("Failure! Could not find any bump in this dataset. Exiting... \n", sep="")
        used <- NULL

        # Cross-validated maximum peeling length from all replicates
        CV.maxsteps <- NULL
        # Cross-validated optimal peeling length over the replicates
        CV.nsteps <- NULL
        # Common selected covariates over the replicates
        CV.sel <- NULL
        # Modal sign values of selected covariates over the replicates (by step)
        CV.sign <- NULL
        # List of CV profiles
        CV.profile <- NULL

  }

  # Creating the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                 "cv.nsteps"=CV.nsteps,
                 "cv.model"=CV.model,
                 "cv.nmodel"=M,
                 "cv.sel"=CV.sel,
                 "cv.sign"=CV.sign,
                 "cv.profile"=CV.profile)
  cat("Cross-Validation and variable selection of Survival Bump Hunting model finished!\n")

  # Returning the final object of class 'CV'
  return(structure(list("x"=x, "times"=times, "status"=status,
                        "B"=B,
                        "K"=K,
                        "cvtype"=cvtype,
                        "cvcriterion"=cvcriterion,
                        "conservative"=conservative,
                        "cvfit"=CV.fit,
                        "arg"=arg,
                        "fdr"=fdr,
                        "thr"=thr,
                        "success"=success,
                        "seed"=seed),
                   class="CV"))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                   sbh(cvobj=NULL,
#                       A=1000, cpv=FALSE, decimals=2,
#                       probval=NULL, timeval=NULL,
#                       parallel=FALSE, conf=NULL, seed=NULL)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

sbh <- function(cvobj=NULL,
                A=1000, cpv=FALSE, decimals=2,
                probval=NULL, timeval=NULL,
                parallel=FALSE, conf=NULL, seed=NULL) {

  # Checking
  if (is.null(cvobj)) {
    stop("Must provide a CV object of the Survival Bump Hunting model. Exiting...\n\n")
  } else {
    if (!inherits(cvobj, 'CV'))
        stop("Primary argument much be an object of class 'CV'")
    if (cvobj$success) {
        # Retrieving data
        x <- cvobj$x
        times <- cvobj$times
        status <- cvobj$status

        # Retrieving CV estimates
        CV.maxsteps <- cvobj$cvfit$cv.maxsteps
        CV.nsteps <- cvobj$cvfit$cv.nsteps
        CV.model <- cvobj$cvfit$cv.model
        CV.sel <- cvobj$cvfit$cv.sel[[CV.model]]
        CV.sign <- cvobj$cvfit$cv.sign[[CV.model]]

        # Pre-selected covariates
        x.sel <- x[, CV.sel, drop=FALSE]
        p.sel <- length(CV.sel)
        n <- nrow(x.sel)

        # Retrieving, parsing and evaluating parameters
        B <- cvobj$B
        K <- cvobj$K
        cvtype <- cvobj$cvtype
        cvcriterion <- cvobj$cvcriterion
        conservative <- cvobj$conservative
        alpha <- NULL
        beta <- NULL
        minn <- NULL
        peelcriterion <- NULL
        eval(parse( text=unlist(strsplit(x=cvobj$arg, split=",")) ))
        arg <- paste("beta=", beta, ",alpha=", alpha, ",minn=", minn, ",L=", CV.nsteps-1, ",peelcriterion=\"", peelcriterion, "\"", sep="")
    } else {
        stop("Submitted CV object previously failed to pre-select any informative covariate or fit the Survival Bump Hunting model. Exiting...\n\n")
    }
  }

  # Review of user choices
  if (is.null(cvobj$fdr) && is.null(cvobj$thr)) {
    cat("No variable selection requested \n")
  } else {
    if (is.null(cvobj$thr)) {
        cat("Requested FDR-based variable selection  \n")
    } else if (is.null(cvobj$fdr)) {
        cat("Requested model-size-based variable selection  \n")
    }
  }

  if ((is.null(cvtype)) || (cvtype == "none") || (is.null(cvcriterion)) || (cvcriterion == "none")) {
    cvtype <- "none"
    cvcriterion <- "none"
    B <- 1
    K <- 1
    conf <- NULL
    parallel <- FALSE
    cat("No cross-validation requested. No replication will be performed. No need of parallelization. \n")
  } else {
    if (B > 1) {
        if (parallel) {
            cat("Requested parallel replicated ", K, "-fold cross-validation with ", conf$cpus*ceiling(B/conf$cpus), " replications \n", sep="")
        } else {
            cat("Requested replicated ", K, "-fold cross-validation with ", B, " replications \n", sep="")
        }
    } else {
        cat("Requested single ", K, "-fold cross-validation without replications \n", sep="")
    }
  }
  cat("Variable pre-selection conservativeness: ", disp(x=conservative), "\n")
  cat("Cross-validation technique: ", disp(x=cvtype), "\n")
  cat("Cross-validation criterion: ", disp(x=cvcriterion), "\n")
  cat("Computation of permutation p-values:", cpv, "\n")
  cat("Peeling criterion: ", disp(x=peelcriterion), "\n")
  cat("Parallelization:", parallel, "\n")
  cat("\n")

  cat("Computation of cross-validated survival estimates from the fitted Survival Bump Hunting model using the PRSP algorithm ... \n")
  if (!parallel) {
    if (is.null(seed)) {
        seed <- runif(n=B, min=1, max=2) * 10^(digits-2)
    } else {
        seed <- (0:(B-1)) + seed
    }
    CV.box.rep.obj <- cv.box.rep(x=x.sel, times=times, status=status,
                                 B=B, K=K, arg=arg,
                                 cvtype=cvtype, decimals=decimals,
                                 varsel=CV.sel, varsign=CV.sign,
                                 probval=probval, timeval=timeval,
                                 parallel=parallel, seed=seed)
  } else {
    if (conf$type == "SOCK") {
        cl <- makeCluster(spec=conf$names,
                          type=conf$type,
                          homogeneous=conf$homo,
                          outfile=conf$outfile,
                          verbose=conf$verbose)
    } else {
        cl <- makeCluster(spec=conf$cpus,
                          type=conf$type,
                          homogeneous=conf$homo,
                          outfile=conf$outfile,
                          verbose=conf$verbose)
    }
    clusterSetRNGStream(cl=cl, iseed=seed)
    a <- ceiling(B/conf$cpus)
    B <- a*conf$cpus
######################################################################################################
clusterEvalQ(cl=cl, expr=library("survival"))
clusterEvalQ(cl=cl, expr=library("Hmisc"))
clusterExport(cl=cl,
              varlist=list("cv.tune.rep",
                           "cv.ave.box", "cv.comb.box",
                           "cv.ave.fold", "cv.comb.fold",
                           "cv.ave.peel", "cv.comb.peel",
                           "cv.folds", "peel.box", "prsp",
                           "updatecut", "endpoints",
                           "is.empty", "cbindlist", "list2mat", "list2array",
                           "lapply.mat", "lapply.array"),
              envir=.GlobalEnv)
######################################################################################################
    obj.cl <- clusterCall(cl=cl, fun=cv.box.rep,
                          x=x.sel, times=times, status=status,
                          B=a, K=K, arg=arg,
                          cvtype=cvtype, decimals=decimals,
                          varsel=CV.sel, varsign=CV.sign,
                          probval=probval, timeval=timeval,
                          parallel=parallel, seed=NULL)
    stopCluster(cl)
    CV.box.rep.obj <- list("cv.trace"=vector(mode="list", length=B),
                           "cv.boxind"=vector(mode="list", length=B),
                           "cv.boxcut"=vector(mode="list", length=B),
                           "cv.support"=vector(mode="list", length=B),
                           "cv.lhr"=vector(mode="list", length=B),
                           "cv.lrt"=vector(mode="list", length=B),
                           "cv.cer"=vector(mode="list", length=B),
                           "cv.time.bar"=vector(mode="list", length=B),
                           "cv.prob.bar"=vector(mode="list", length=B),
                           "cv.max.time.bar"=vector(mode="list", length=B),
                           "cv.min.prob.bar"=vector(mode="list", length=B))
    for (b in 1:conf$cpus) {
        CV.box.rep.obj$cv.trace[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.trace
        CV.box.rep.obj$cv.boxind[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.boxind
        CV.box.rep.obj$cv.boxcut[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.boxcut
        CV.box.rep.obj$cv.support[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.support
        CV.box.rep.obj$cv.lhr[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.lhr
        CV.box.rep.obj$cv.lrt[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.lrt
        CV.box.rep.obj$cv.cer[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.cer
        CV.box.rep.obj$cv.time.bar[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.time.bar
        CV.box.rep.obj$cv.prob.bar[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.prob.bar
        CV.box.rep.obj$cv.max.time.bar[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.max.time.bar
        CV.box.rep.obj$cv.min.prob.bar[((b-1)*a+1):(b*a)] <- obj.cl[[b]]$cv.min.prob.bar
    }
    CV.box.rep.obj$success <- obj.cl[[1]]$success
  }
  success <- CV.box.rep.obj$success

  if (success) {

    # Collect the peeling statistics for each step from all the replicates
    CV.trace <- CV.box.rep.obj$cv.trace
    CV.boxind <- CV.box.rep.obj$cv.boxind
    CV.boxcut <- CV.box.rep.obj$cv.boxcut
    CV.support <- CV.box.rep.obj$cv.support
    CV.lhr <- CV.box.rep.obj$cv.lhr
    CV.lrt <- CV.box.rep.obj$cv.lrt
    CV.cer <- CV.box.rep.obj$cv.cer
    CV.time.bar <- CV.box.rep.obj$cv.time.bar
    CV.prob.bar <- CV.box.rep.obj$cv.prob.bar
    CV.max.time.bar <- CV.box.rep.obj$cv.max.time.bar
    CV.min.prob.bar <- CV.box.rep.obj$cv.min.prob.bar

    # Box membership indicator vector of all observations for each step using the modal or majority vote value over the replicates
    cat("Generating cross-validated box memberships for each step ...\n")
    CV.boxind <- lapply.array(X=CV.boxind, rowtrunc=CV.nsteps, FUN=function(x){mean(x, na.rm=TRUE) >= 0.5}, MARGIN=1:2)
    rownames(CV.boxind) <- paste("step", 0:(CV.nsteps-1), sep="")
    colnames(CV.boxind) <- rownames(x.sel)

    # Modal trace values (over the replicates) of covariate usage at each step
    cat("Generating cross-validated modal trace values of covariate usage at each step ...\n")
    trace.dist <- lapply.array(X=CV.trace,
                               rowtrunc=CV.nsteps,
                               FUN=function(x){if (all(is.na(x)))
                                                 return(NA)
                                               else
                                                 return(as.numeric(names(which.max(table(x, useNA="no")))))
                                              },
                               MARGIN=c(1,3))
    dimnames(trace.dist) <- list(paste("step", 0:(CV.nsteps-1), sep=""), 1:B)
    trace.mode <- apply(X=trace.dist,
                        FUN=function(x){as.numeric(names(which.max(table(x, useNA="no"))))},
                        MARGIN=1)
    mat <- pmatch(x=names(CV.sel)[trace.mode], table=colnames(x), nomatch=NA, duplicates.ok=TRUE)
    CV.trace <- c(0, mat[!is.na(mat)])
    names(CV.trace) <- paste("step", 0:(CV.nsteps-1), sep="")

    # Used covariates for peeling at each step, based on covariate trace modal values
    CV.used <- sort(unique(as.numeric(CV.trace[-1])))
    names(CV.used) <- colnames(x)[CV.used]
    cat("Covariates used for peeling at each step, based on covariate trace modal values:\n")
    print(CV.used)

    # Box rules for the pre-selected covariates at each step
    cat("Generating cross-validated box rules for the pre-selected covariates at each step ...\n")
    CV.boxcut.mu <- round(lapply.array(X=CV.boxcut, rowtrunc=CV.nsteps, FUN=function(x){mean(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
    CV.boxcut.sd <- round(lapply.array(X=CV.boxcut, rowtrunc=CV.nsteps, FUN=function(x){sd(x, na.rm=TRUE)}, MARGIN=1:2), digits=decimals)
    rownames(CV.boxcut.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
    rownames(CV.boxcut.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
    colnames(CV.boxcut.mu) <- names(CV.sel)
    colnames(CV.boxcut.sd) <- names(CV.sel)
    CV.frame <- as.data.frame(matrix(data=NA, nrow=CV.nsteps, ncol=p.sel, dimnames=list(paste("step", 0:(CV.nsteps-1), sep=""),names(CV.sel))))
    for (j in 1:p.sel) {
        if (CV.sign[j] > 0) {
            ss <- ">="
        } else {
            ss <- "<="
        }
        CV.frame[, j] <- paste(paste(names(CV.sel)[j], ss, format(x=CV.boxcut.mu[, j], digits=decimals, nsmall=decimals), sep=""),
                               format(x=CV.boxcut.sd[, j], digits=decimals, nsmall=decimals), sep=" +/- ")
    }
    CV.rules <- list("mean"=CV.boxcut.mu, "sd"=CV.boxcut.sd, "frame"=CV.frame)

    # Box statistics for each step
    cat("Generating cross-validated box statistics for each step ...\n")
    CV.support.mu <- round(lapply.mat(X=CV.support, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.support.sd <- round(lapply.mat(X=CV.support, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.size.mu <- round(n*CV.support.mu,0)
    CV.size.sd <- n*CV.support.sd
    CV.lhr.mu <- round(lapply.mat(X=CV.lhr, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.lhr.sd <- round(lapply.mat(X=CV.lhr, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.lrt.mu <- round(lapply.mat(X=CV.lrt, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.lrt.sd <- round(lapply.mat(X=CV.lrt, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.cer.mu <- round(lapply.mat(X=CV.cer, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.cer.sd <- round(lapply.mat(X=CV.cer, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.time.bar.mu <- round(lapply.mat(X=CV.time.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.time.bar.sd <- round(lapply.mat(X=CV.time.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.prob.bar.mu <- round(lapply.mat(X=CV.prob.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.prob.bar.sd <- round(lapply.mat(X=CV.prob.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.max.time.bar.mu <- round(lapply.mat(X=CV.max.time.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.max.time.bar.sd <- round(lapply.mat(X=CV.max.time.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.min.prob.bar.mu <- round(lapply.mat(X=CV.min.prob.bar, FUN=function(x){mean(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.min.prob.bar.sd <- round(lapply.mat(X=CV.min.prob.bar, FUN=function(x){sd(x, na.rm=TRUE)}, coltrunc=CV.nsteps), digits=decimals)
    CV.stats.mu <- data.frame("cv.support"=CV.support.mu,
                              "cv.size"=CV.size.mu,
                              "cv.lhr"=CV.lhr.mu,
                              "cv.lrt"=CV.lrt.mu,
                              "cv.cer"=CV.cer.mu,
                              "cv.time.bar"=CV.time.bar.mu,
                              "cv.prob.bar"=CV.prob.bar.mu,
                              "cv.max.time.bar"=CV.max.time.bar.mu,
                              "cv.min.prob.bar"=CV.min.prob.bar.mu)
    rownames(CV.stats.mu) <- paste("step", 0:(CV.nsteps-1), sep="")
    CV.stats.sd <- data.frame("cv.support"=CV.support.sd,
                              "cv.size"=CV.size.sd,
                              "cv.lhr"=CV.lhr.sd,
                              "cv.lrt"=CV.lrt.sd,
                              "cv.cer"=CV.cer.sd,
                              "cv.time.bar"=CV.time.bar.sd,
                              "cv.prob.bar"=CV.prob.bar.sd,
                              "cv.max.time.bar"=CV.max.time.bar.sd,
                              "cv.min.prob.bar"=CV.min.prob.bar.sd)
    rownames(CV.stats.sd) <- paste("step", 0:(CV.nsteps-1), sep="")
    CV.stats <- list("mean"=CV.stats.mu, "sd"=CV.stats.sd)

    # P-values for each step
    if (cpv) {
        cat("Computation of cross-validated permutation p-values for each step ... \n")
        arg <- paste("beta=", beta, ",alpha=", alpha, ",minn=", minn, ",L=", CV.nsteps-1, ",peelcriterion=\"", peelcriterion, "\"", sep="")
        CV.pval <- cv.pval(x=x.sel, times=times, status=status,
                           cvtype=cvtype, decimals=decimals,
                           varsel=CV.sel, varsign=CV.sign,
                           A=A, K=K, arg=arg,
                           obs.chisq=CV.stats$mean$cv.lrt,
                           parallel=parallel, conf=conf)
    } else {
        CV.pval <- NULL
    }

  } else {

    # Selected covariates
    cat("Failed to pre-select any informative covariate after 10 successive trials. Exiting...\n", sep="")
    # Cross-validated maximum peeling length from all replicates
    CV.maxsteps <- NULL
    # Cross-validated optimal length from all replicates
    CV.nsteps <- NULL
    # Modal or majority vote trace value over the replicates
    CV.trace <- NULL
    # Common selected covariates over the replicates
    CV.sel <- NULL
    # Trace values of selected covariates over the replicates
    CV.sign <- NULL
    # Box boundaries and box peeling rules for each step
    CV.rules <- NULL
    # Box membership indicator vector of all observations for each step
    CV.boxind <- NULL
    # Box statistics for each step
    CV.stats <- NULL
    # P-values for each step
    CV.pval <- NULL

  }

  # Creating the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.maxsteps,
                 "cv.nsteps"=CV.nsteps,
                 "cv.trace"=CV.trace,
                 "cv.sign"=CV.sign,
                 "cv.sel"=CV.sel,
                 "cv.used"=CV.used,
                 "cv.boxind"=CV.boxind,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.pval"=CV.pval)
  cat("Computation of survival estimates finished!\n")

  # Returning the final object of class 'PRSP'
  return(structure(list("x"=x, "times"=times, "status"=status,
                        "B"=B, "K"=K, "A"=A, "cpv"=cpv,
                        "decimals"=decimals, "arg"=arg,
                        "cvtype"=cvtype,
                        "cvcriterion"=cvcriterion,
                        "conservative"=conservative,
                        "probval"=probval, "timeval"=timeval,
                        "cvfit"=CV.fit,
                        "plot"=success,
                        "config"=list("parallel"=parallel,
                                      "names"=conf$names,
                                      "cpus"=conf$cpus,
                                      "type"=conf$type,
                                      "homo"=conf$homo,
                                      "verbose"=conf$verbose,
                                      "outfile"=conf$outfile),
                        "seed"=seed),
                   class="PRSP"))
}
##########################################################################################################################################




##########################################################################################################################################
# 2. END-USER FUNCTION FOR NEWS
##########################################################################################################################################

##########################################################################################################################################
################
#Usage         :
################
#                   PRIMsrc.news(...)
#
################
# Description   :
################
#                   Function to display the log file of updates of the PRIMsrc package.
#
################
# Arguments     :
################
# ...               Further arguments passed to or from other methods.
#
################
# Values        :
################
#
##########################################################################################################################################

PRIMsrc.news <- function(...) {
    newsfile <- file.path(system.file(package="PRIMsrc"), "NEWS")
    file.show(newsfile)
}
##########################################################################################################################################




##########################################################################################################################################
# 3. END-USER S3-GENERIC FUNCTIONS FOR SUMMARY, PRINT, PLOT AND PREDICTION
##########################################################################################################################################

##########################################################################################################################################
################
#Usage         :
################
#                   summary(object, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

summary.PRSP <- function(object, ...) {

  if (!inherits(object, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  alpha <- NULL
  beta <- NULL
  minn <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=object$arg, split=",")) ))

  cat("S3-class object: `", attr(x=object, "class"), "` \n\n")

  if (object$cvtype != "none") {
    if (object$B > 1) {
        if (object$config$parallel) {
            cat("Replicated ", object$K, "-fold cross-validation with ", object$config$cpus*ceiling(object$B/object$config$cpus), " replications \n\n", sep="")
        } else {
            cat("Replicated ", object$K, "-fold cross-validation with ", object$B, " replications \n\n", sep="")
        }
    } else {
        cat("Single ", object$K, "-fold cross-validation without replications \n\n", sep="")
    }
  } else {
    cat("'PRSP' object without cross-validation and replications\n\n", sep="")
  }
  cat("Variable pre-selection conservativeness: ", disp(x=object$conservative), "\n\n")
  cat("PRSP parameters:\n")
  cat("\t Peeling criterion: ", disp(x=peelcriterion), "\n")
  cat("\t Peeling percentile: ", alpha*100, "%\n")
  cat("\t Minimal box support: ", beta*100, "%\n")
  cat("\t Minimal box sample size: ", minn, "\n\n")
  cat("Cross-validation technique: ", disp(x=object$cvtype), "\n\n")
  cat("Cross-validation criterion: ", disp(x=object$cvcriterion), "\n\n")
  cat("Number of decimals: ", object$decimals, "\n\n")
  cat("Computation of permutation p-values: ", object$cpv, "\n\n")
  cat("Parallelization: ", object$config$parallel, "\n\n")

  invisible()
}
##########################################################################################################################################




##########################################################################################################################################
################
#Usage         :
################
#                   print(x, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

print.PRSP <- function(x, ...) {

  if (!inherits(x, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  obj <- x
  cat("Selected covariates:\n")
  print(obj$cvfit$cv.sel)
  cat("\n")

  cat("Used covariates:\n")
  print(obj$cvfit$cv.used)
  cat("\n")

  cat("Maximum number of peeling steps (counting step #0):\n")
  print(obj$cvfit$cv.maxsteps)
  cat("\n")

  out <- obj$cvfit$cv.nsteps
  names(out) <- NULL
  cat("Optimal number of peeling steps (counting step #0):\n")
  print(out)
  cat("\n")

  cat("Modal trace values of covariate usage at each peeling step:\n")
  print(obj$cvfit$cv.trace)
  cat("\n")

  cat("Cross-validated permutation p-values at each peeling step:\n")
  print(obj$cvfit$cv.pval, quote = FALSE)
  cat("\n")

  cat("Decision rules for the used covariates (columns) at each peeling step (rows):\n")
  print(obj$cvfit$cv.rules$frame[,names(obj$cvfit$cv.used),drop=FALSE], quote = FALSE)
  cat("\n")

  colnames(obj$cvfit$cv.stats$mean) <- c("Support", "Size", "LHR", "LRT", "CER", "EFT", "EFP", "MEFT", "MEFP")
  cat("Box endpoint quantities of interest (columns) at each peeling step (rows):\n")
  print(obj$cvfit$cv.stats$mean)
  cat("\n")

  cat("Individual observation box membership indicator (columns) at each peeling step (rows):\n")
  print(obj$cvfit$cv.boxind)
  cat("\n")

  invisible()
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    plot(x,
#                         main=NULL,
#                         proj=c(1,2), splom=TRUE, boxes=FALSE,
#                         steps=x$cvfit$cv.nsteps,
#                         pch=16, cex=0.5, col=, col=2:(length(steps)+1),
#                         col.box=2:(length(steps)+1), lty.box=rep(2,length(steps)), lwd.box=rep(1,length(steps)),
#                         add.legend=TRUE,
#                         device=NULL, file="Scatter Plot", path=getwd(),
#                         horizontal=FALSE, width=5, height=5, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot.PRSP <- function(x,
                      main=NULL,
                      proj=c(1,2), splom=TRUE, boxes=FALSE,
                      steps=x$cvfit$cv.nsteps,
                      pch=16, cex=0.5, col=2:(length(steps)+1),
                      col.box=2:(length(steps)+1), lty.box=rep(2,length(steps)), lwd.box=rep(1,length(steps)),
                      add.legend=TRUE,
                      device=NULL, file="Scatter Plot", path=getwd(),
                      horizontal=FALSE, width=5, height=5, ...) {

    if (!inherits(x, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

    obj <- x
    if (obj$plot) {

        scatterplot <- function(obj,
                                main,
                                proj, splom, boxes,
                                steps,
                                add.legend, pch, cex, col,
                                col.box, lty.box, lwd.box, ...) {

        if (!is.null(main)) {
            par(mfrow=c(1, 1), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
        } else {
            par(mfrow=c(1, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
        }

        toplot <- obj$cvfit$cv.used[proj]
        varnames <- colnames(obj$x)
        X <- obj$x[,varnames[toplot],drop=FALSE]
        X.names <- colnames(X)

        if (is.null(steps))
          steps <- obj$cvfit$cv.nsteps

        L <- length(steps)
        eqscplot(x=X, type="p", pch=pch, cex=cex, col=1, main=NULL, xlab=X.names[1], ylab=X.names[2], ...)
        if (splom) {
            for (i in 1:L) {
                w <- obj$cvfit$cv.boxind[steps[i],]
                points(x=obj$x[w,varnames[toplot],drop=FALSE], type="p", pch=pch, cex=cex, col=col[i], ...)
            }
        }
        if (boxes) {
            X.range <- apply(X=X, MARGIN=2, FUN=range)
            boxcut <- obj$cvfit$cv.rules$mean[steps,varnames[toplot],drop=FALSE]
            varsign <- obj$cvfit$cv.sign[varnames[toplot]]
            vertices <- vector(mode="list", length=L)
            for (i in 1:L) {
                vertices[[i]] <- matrix(data=NA, nrow=2, ncol=2, dimnames=list(c("LB","UB"), X.names))
                for (j in 1:2) {
                    vertices[[i]][1,j] <- ifelse(test=(varsign[j] > 0),
                                                 yes=max(X.range[1,j], boxcut[i,j]),
                                                 no=min(X.range[1,j], boxcut[i,j]))
                    vertices[[i]][2,j] <- ifelse(test=(varsign[j] < 0),
                                                 yes=min(X.range[2,j], boxcut[i,j]),
                                                 no=max(X.range[2,j], boxcut[i,j]))
                }
            }
            for (i in 1:L) {
                rect(vertices[[i]][1,1], vertices[[i]][1,2], vertices[[i]][2,1], vertices[[i]][2,2],
                     border=col.box[i], col=NA, lty=lty.box[i], lwd=lwd.box[i])
            }
        }
        if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
        }
        if (add.legend) {
            legend("topleft", xpd=TRUE, inset=0.01, legend=paste("Step: ", steps, sep=""), pch=pch, col=col, cex=cex)
        }
    }

    if (is.null(device)) {
        scatterplot(obj=obj,
                    main=main,
                    proj=proj, splom=splom, boxes=boxes, steps=steps,
                    add.legend=add.legend, pch=pch, cex=cex, col=col,
                    col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        scatterplot(obj=obj,
                    main=main,
                    proj=proj, splom=splom, boxes=boxes, steps=steps,
                    add.legend=add.legend, pch=pch, cex=cex, col=col,
                    col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        scatterplot(obj=obj,
                    main=main,
                    proj=proj, splom=splom, boxes=boxes, steps=steps,
                    add.legend=add.legend, pch=pch, cex=cex, col=col,
                    col.box=col.box, lty.box=lty.box, lwd.box=lwd.box)
        dev.off()
    } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
  } else {
    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")
  }
  invisible()
}
##########################################################################################################################################




##########################################################################################################################################
################
#Usage         :
################
#                   predict(object, newdata, steps, na.action = na.omit, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

predict.PRSP <- function (object, newdata, steps, na.action = na.omit, ...) {

  if (!inherits(object, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  X <- as.matrix(newdata)
  X.names <- colnames(X)
  X.range <- apply(X=X, MARGIN=2, FUN=range)
  n <- nrow(X)
  p <- ncol(X)

  toplot <- object$cvfit$cv.used
  varnames <- colnames(object$x)

  if (length(toplot) != p) {
    stop("Non-matching dimensionality of newdata to PRSP object of used covariates.\n")
  }

  if (missing(steps) || is.null(steps))
    steps <- object$cvfit$cv.nsteps

  L <- length(steps)
  boxcut <- object$cvfit$cv.rules$mean[steps,varnames[toplot],drop=FALSE]
  varsign <- object$cvfit$cv.sign[varnames[toplot]]

  pred.boxind <- matrix(NA, nrow=L, ncol=n, dimnames=list(paste("step ", steps, sep=""), rownames(X)))
  for (l in 1:L) {
    boxcutsign <- boxcut[l, ] * varsign
    x.cut <- t(t(X) * varsign)
    x.ind <- t(t(x.cut) >= boxcutsign)
    pred.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boundaries for all axes directions
  }

  pred.vertices <- vector(mode="list", length=L)
  names(pred.vertices) <- paste("step ", steps, sep="")
  for (i in 1:L) {
    pred.vertices[[i]] <- matrix(data=NA, nrow=2, ncol=p, dimnames=list(c("LB","UB"), X.names))
    for (j in 1:p) {
      pred.vertices[[i]][1,j] <- ifelse(test=(varsign[j] > 0),
                                        yes=max(X.range[1,j], boxcut[i,j]),
                                        no=min(X.range[1,j], boxcut[i,j]))
      pred.vertices[[i]][2,j] <- ifelse(test=(varsign[j] < 0),
                                        yes=min(X.range[2,j], boxcut[i,j]),
                                        no=max(X.range[2,j], boxcut[i,j]))
    }
  }

  return(list("boxind"=pred.boxind, "vertices"=pred.vertices))
}
##########################################################################################################################################




##########################################################################################################################################
# 4. END-USER PLOTTING FUNCTIONS FOR MODEL VALIDATION AND VISUALIZATION OF RESULTS
##########################################################################################################################################

##########################################################################################################################################
################
# Usage         :
################
#                    plot_profile(cvobj,
#                                 main=NULL,
#                                 xlab="Peeling Steps", ylab="Mean Profiles",
#                                 add.sd=TRUE, add.profiles=TRUE,
#                                 pch=20, col=1, lty=1, lwd=2, cex=2,
#                                 device=NULL, file="Profile Plot", path=getwd(),
#                                 horizontal=FALSE, width=8.5, height=11.5, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot_profile <- function(cvobj,
                         main=NULL,
                         xlab="Peeling Steps", ylab="Mean Profiles",
                         add.sd=TRUE, add.profiles=TRUE,
                         pch=20, col=1, lty=1, lwd=2, cex=2,
                         device=NULL, file="Profile Plot", path=getwd(),
                         horizontal=FALSE, width=8.5, height=11.5, ...) {

  if (!inherits(cvobj, 'CV'))
        stop("Primary argument much be a CV object")

  if (cvobj$success) {
    if (cvobj$cvtype == "none") {

      cat("No CV here, so no cross-validated tuning profile to plot!\n")

    } else {

      profileplot <- function(cvobj, main, xlab, ylab,
                              add.sd, add.profiles,
                              pch, col, lty, lwd, cex, ...) {

        if (is.null(cvobj$fdr) && is.null(cvobj$thr)) {
            model.size <- ncol(cvobj$x)
        } else {
            if (is.null(cvobj$thr)) {
                model.size <- cvobj$fdr
            } else if (is.null(cvobj$fdr)) {
                model.size <- cvobj$thr
            }
        }
        M <- length(model.size)
        profiles <- cvobj$cvfit$cv.profile$profiles
        mu.profile <- cvobj$cvfit$cv.profile$mean
        se.profile <- cvobj$cvfit$cv.profile$se
        if (cvobj$cvcriterion == "lhr") {
          txt <- "LHR"
          optisteps <- apply(X=mu.profile, MARGIN=1, which.max)
        } else if (cvobj$cvcriterion == "lrt") {
          txt <- "LRT"
          optisteps <- apply(X=mu.profile, MARGIN=1, which.max)
        } else if (cvobj$cvcriterion == "cer") {
          txt <- "CER"
          optisteps <- apply(X=mu.profile, MARGIN=1, which.min)
        } else {
          stop("Invalid CV criterion.\n")
        }
        subplots <- min(M,5)
        if (!is.null(main)) {
          par(mfrow=c(subplots, 1), oma=c(0, 0, 2, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
        } else {
          par(mfrow=c(subplots, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 4.0, 1.5), mgp=c(1.5, 0.5, 0))
        }
        Lm <- cvobj$cvfit$cv.maxsteps
        for (m in 1:M) {
            ylim <- range(0, 1, mu.profile[m,] - se.profile[m,], mu.profile[m,] + se.profile[m,], na.rm=TRUE)
            if (add.profiles) {
                matplot(profiles[m,,], axes=FALSE, type="b",
                        xlab="", ylab="",
                        main="", ylim=ylim,
                        pch=pch, lty=1, lwd=lwd/4, cex=cex/4, ...)
                par(new=TRUE)
            }
            plot(0:(Lm-1), mu.profile[m,], axes=FALSE, type="b",
                 xlab=xlab, ylab=paste(txt ," ", ylab, sep=""),
                 main=paste("Model size: ", model.size[m], sep=""), ylim=ylim,
                 pch=pch, col=col, lty=lty, lwd=lwd, cex=cex, ...)
            axis(side=1, pos=min(ylim), at=0:(Lm-1), labels=0:(Lm-1), cex.axis=1, line=NA)
            axis(side=2, pos=0, at=pretty(ylim), cex.axis=1, line=NA)
            segments(x0=optisteps[m]-1, y0=min(ylim), x1=optisteps[m]-1, y1=mu.profile[m,optisteps[m]], col=col, lty=2, lwd=lwd)
            if (add.sd) {
                arrows(0:(Lm-1), mu.profile[m,], 0:(Lm-1), mu.profile[m,] - se.profile[m,], length=0.1, angle=90, code=2, col=col, lwd=lwd/2)
                arrows(0:(Lm-1), mu.profile[m,], 0:(Lm-1), mu.profile[m,] + se.profile[m,], length=0.1, angle=90, code=2, col=col, lwd=lwd/2)
            }
        }
        if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, line=1, outer=TRUE)
        }
      }
    }

    if (is.null(device)) {
        cat("Device: ",  dev.cur(), "\n")
        profileplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    add.sd=add.sd, add.profiles=add.profiles,
                    pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        cat("Device: ",  dev.cur(), "\n")
        profileplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    add.sd=add.sd, add.profiles=add.profiles,
                    pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        cat("Device: ",  dev.cur(), "\n")
        profileplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    add.sd=add.sd, add.profiles=add.profiles,
                    pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
        dev.off()
    } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }

  } else {

    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")

  }
  invisible()
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                    plot_surface(cvobj,
#                                 main=NULL,
#                                 xlab="Model Size (# variables)", ylab="Peeling Steps",
#                                 theta=5, phi=10, expand=0.2, col.surf="lightblue",
#                                 add.line=FALSE, col="yellow", pch=20, lty=1, lwd=1, cex=1,
#                                 device=NULL, file="Surface Plot", path=getwd(),
#                                 horizontal=FALSE, width=8.5, height=5.0, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot_surface <- function(cvobj,
                         main=NULL,
                         xlab="Model Size (# variables)", ylab="Peeling Steps",
                         theta=5, phi=10, expand=0.2, col.surf="lightblue",
                         add.line=FALSE, col="yellow", pch=20, lty=1, lwd=1, cex=1,
                         device=NULL, file="Surface Plot", path=getwd(),
                         horizontal=FALSE, width=8.5, height=5.0, ...) {

  if (!inherits(cvobj, 'CV'))
        stop("Primary argument much be a CV object")

  if (cvobj$success) {
    if (cvobj$cvtype == "none") {

      cat("No CV here, so no cross-validated tuning profile to plot!\n")

    } else if (is.null(cvobj$fdr) && is.null(cvobj$thr)) {

      cat("No model selection here, so no mean CV surface to plot!\n")

    } else {

       surfaceplot <- function(cvobj, main, xlab, ylab,
                               theta, phi, expand, col.surf,
                               add.line, pch, col, lty, lwd, cex, ...) {

        if (is.null(cvobj$thr)) {
            models <- cvobj$fdr
        } else if (is.null(cvobj$fdr)) {
            models <- cvobj$thr
        }
        M <- length(models)
        optima <- numeric(M)
        optisteps <- numeric(M)
        steps <- seq(cvobj$cvfit$cv.maxsteps)
        mu.profile <- cvobj$cvfit$cv.profile$mean
        se.profile <- cvobj$cvfit$cv.profile$se

        if (cvobj$cvcriterion == "lhr") {
          txt <- "LHR"
          for (m in 1:M) {
            optisteps[m] <- which.max(mu.profile[m,])
            optima[m] <- mu.profile[m,optisteps[m]]
          }
        } else if (cvobj$cvcriterion == "lrt") {
          txt <- "LRT"
          for (m in 1:M) {
            optisteps[m] <- which.max(mu.profile[m,])
            optima[m] <- mu.profile[m,optisteps[m]]
          }
        } else if (cvobj$cvcriterion == "cer") {
          txt <- "CER"
          for (m in 1:M) {
            optisteps[m] <- which.min(mu.profile[m,])
            optima[m] <- mu.profile[m,optisteps[m]]
          }
        } else {
          stop("Invalid CV criterion.\n")
        }
        res <- persp(x=models, y=steps, z=mu.profile,
                     theta=theta, phi=phi, expand=expand,
                     col=col.surf, shade=TRUE,
                     xlab=xlab, ylab=ylab, zlab=txt, main=main,
                     ticktype = "detailed", box=TRUE, nticks=10, ...)
        if (add.line) {
            points(trans3d(x=models, y=optisteps, z=optima, pmat=res), col=col, pch=pch, cex=cex, ...)
            lines(trans3d(x=models, y=optisteps, z=optima, pmat=res), col=col, lwd=lwd, lty=lty, ...)
        }
        if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
        }
      }

      if (is.null(device)) {
        cat("Device: ",  dev.cur(), "\n")
        surfaceplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    theta=theta, phi=phi, expand=expand, col.surf=col.surf,
                    add.line=add.line, pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
      } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        cat("Device: ",  dev.cur(), "\n")
        surfaceplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    theta=theta, phi=phi, expand=expand, col.surf=col.surf,
                    add.line=add.line, pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
        dev.off()
      } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        cat("Device: ",  dev.cur(), "\n")
        surfaceplot(cvobj=cvobj, main=main, xlab=xlab, ylab=ylab,
                    theta=theta, phi=phi, expand=expand, col.surf=col.surf,
                    add.line=add.line, pch=pch, col=col, lty=lty, lwd=lwd, cex=cex)
        dev.off()
      } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
      }
    }

  } else {

    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")

  }
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    plot_boxtraj (object,
#                                  main=NULL,
#                                  xlab="Box Mass", ylab="Covariate Range",
#                                  toplot=object$cvfit$cv.used,
#                                  col.cov, lty.cov, lwd.cov,
#                                  col=1, lty=1, lwd=1,
#                                  cex=1, add.legend=FALSE, text.legend=NULL,
#                                  nr=NULL, nc=NULL,
#                                  device=NULL, file="Covariate Trajectory Plots", path=getwd())
#                                  horizontal=FALSE, width=8.5, height=8.5, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot_boxtraj <- function(object,
                         main=NULL,
                         xlab="Box Mass", ylab="Covariate Range",
                         toplot=object$cvfit$cv.used,
                         col.cov, lty.cov, lwd.cov,
                         col=1, lty=1, lwd=1,
                         cex=1, add.legend=FALSE, text.legend=NULL,
                         nr=NULL, nc=NULL,
                         device=NULL, file="Covariate Trajectory Plots", path=getwd(),
                         horizontal=FALSE, width=8.5, height=11.5, ...) {

  if (!inherits(object, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  if (object$plot) {
    boxtrajplot <- function(object,
                            main, xlab, ylab,
                            toplot,
                            col.cov, lty.cov, lwd.cov,
                            col, lty, lwd,
                            cex, add.legend, text.legend,
                            nr, nc, ...) {
        p <- length(toplot)
        varnames <- colnames(object$x)
        if (is.null(nc))
            nc <- 3
        if (is.null(nr)) {
            if (p %% nc == 0) {
                nr <- p%/%nc + 2
            } else {
                nr <- ((p+(1:nc))[which((p+(1:nc)) %% nc == 0)])%/%nc + 2
            }
        }
        if (missing(col.cov)) {
            col.cov <- 2:(p+1)
        }
        if (missing(lty.cov)) {
            lty.cov <- rep(1,p)
        }
        if (missing(lwd.cov)) {
            lwd.cov <- rep(1,p)
        }
        if (!is.null(main)) {
            par(mfrow=c(nr, nc), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 2.0, 1.5), mgp=c(1.5, 0.5, 0))
        } else {
            par(mfrow=c(nr, nc), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 2.0, 1.5), mgp=c(1.5, 0.5, 0))
        }
        for (j in 1:p) {
            plot(x=object$cvfit$cv.stats$mean$cv.support,
                 y=object$cvfit$cv.rules$mean[,varnames[toplot[j]]],
                 type='s', col=col.cov[j], lty=lty.cov[j], lwd=lwd.cov[j],
                 main=paste(varnames[toplot[j]], " covariate trajectory", sep=""), cex.main=cex,
                 xlim=range(0,1),
                 ylim=range(object$x[,toplot[j]], na.rm=TRUE),
                 xlab=xlab,
                 ylab=ylab, ...)
        }
        if (add.legend)
          legend("bottomleft", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr-1, 1))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.support,
             type='s', col=col, lty=lty, lwd=lwd,
             main="Box support trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0,1),
             xlab=xlab,
             ylab=expression(paste("Support (", beta, ")", sep="")), ...)
        if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr-1, 2))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.max.time.bar,
             type='s', col=col, lty=lty, lwd=lwd,
             main="MEFT trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0, object$cvfit$cv.stats$mean$cv.max.time.bar, na.rm=TRUE),
             xlab=xlab,
             ylab="Time", ...)
        if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr-1, 3))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.min.prob.bar,
             type='s', col=col, lty=lty, lwd=lwd,
             main="MEFP trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0,1),
             xlab=xlab,
             ylab="Probability", ...)
        if (add.legend)
            legend("bottomright", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr, 1))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.lhr,
             type='s', col=col, lty=lty, lwd=lwd,
             main="LHR trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0, object$cvfit$cv.stats$mean$cv.lhr, na.rm=TRUE),
             xlab=xlab,
             ylab=expression(paste("Log-Hazard Ratio (", lambda,")", sep="")), ...)
        if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr, 2))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.lrt,
             type='s', col=col, lty=lty, lwd=lwd,
             main="LRT trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0, object$cvfit$cv.stats$mean$cv.lrt, na.rm=TRUE),
             xlab=xlab,
             ylab=expression(paste("Log-rank test (", chi^2 ,")", sep="")), ...)
        if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
        par(mfg=c(nr, 3))
        plot(object$cvfit$cv.stats$mean$cv.support,
             object$cvfit$cv.stats$mean$cv.cer,
             type='s', col=col, lty=lty, lwd=lwd,
             main="CER trajectory", cex.main=cex,
             xlim=range(0,1),
             ylim=range(0,1),
             xlab=xlab,
             ylab=expression(paste("1-C (", theta,")", sep="")), ...)
        if (add.legend)
            legend("top", inset=0.01, legend=text.legend, cex=cex)
        if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
        }
    }

    if (is.null(device)) {
        cat("Device: ",  dev.cur(), "\n")
        boxtrajplot(object=object,
                    main=main, xlab=xlab, ylab=ylab,
                    toplot=toplot,
                    col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                    col=col, lty=lty, lwd=lwd,
                    cex=cex, add.legend=add.legend, text.legend=text.legend,
                    nr=nr, nc=nc)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        cat("Device: ",  dev.cur(), "\n")
        boxtrajplot(object=object,
                    main=main, xlab=xlab, ylab=ylab,
                    toplot=toplot,
                    col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                    col=col, lty=lty, lwd=lwd,
                    cex=cex, add.legend=add.legend, text.legend=text.legend,
                    nr=nr, nc=nc)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        cat("Device: ",  dev.cur(), "\n")
        boxtrajplot(object=object,
                    main=main, xlab=xlab, ylab=ylab,
                    toplot=toplot,
                    col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                    col=col, lty=lty, lwd=lwd,
                    cex=cex, add.legend=add.legend, text.legend=text.legend,
                    nr=nr, nc=nc)
        dev.off()
    } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
  } else {
    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")
  }
  invisible()
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    plot_boxtrace (object,
#                                   main=NULL,
#                                   xlab="Box Mass", ylab="Covariate Range (centered)",
#                                   toplot=object$cvfit$cv.used,
#                                   center=TRUE, scale=FALSE,
#                                   col.cov, lty.cov, lwd.cov,
#                                   col=1, lty=1, lwd=1,
#                                   cex=1, add.legend=FALSE, text.legend=NULL,
#                                   device=NULL, file="Covariate Trace Plots", path=getwd(),
#                                   horizontal=FALSE, width=8.5, height=8.5, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot_boxtrace <- function(object,
                          main=NULL, xlab="Box Mass", ylab="Covariate Range (centered)",
                          toplot=object$cvfit$cv.used,
                          center=TRUE, scale=FALSE,
                          col.cov, lty.cov, lwd.cov,
                          col=1, lty=1, lwd=1,
                          cex=1, add.legend=FALSE, text.legend=NULL,
                          device=NULL, file="Covariate Trace Plots", path=getwd(),
                          horizontal=FALSE, width=8.5, height=8.5, ...) {

  if (!inherits(object, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  if (object$plot) {
    boxtraceplot <- function(object,
                             main, xlab, ylab,
                             toplot,
                             center, scale,
                             col.cov, lty.cov, lwd.cov,
                             col, lty, lwd,
                             cex, add.legend, text.legend, ...) {
        p <- length(toplot)
        varnames <- colnames(object$x)
        maxlength <- max(sapply(X=varnames, FUN=function(x){nchar(x, type="chars", allowNA=TRUE)}))

        if (missing(col.cov)) {
            col.cov <- 2:(p+1)
        }
        if (missing(lty.cov)) {
            lty.cov <- rep(1,p)
        }
        if (missing(lwd.cov)) {
            lwd.cov <- rep(1,p)
        }
        if (!is.null(main)) {
            par(mfrow=c(2, 1), oma=c(0, 0, 2, 0), mar=c(2.5, 2+maxlength/2, 2.0, 0.0), mgp=c(1.5, 0.5, 0))
        } else {
            par(mfrow=c(2, 1), oma=c(0, 0, 0, 0), mar=c(2.5, 2+maxlength/2, 2.0, 0.0), mgp=c(1.5, 0.5, 0))
        }

        boxcut.scaled <- scale(x=object$cvfit$cv.rules$mean[,varnames[toplot],drop=FALSE], center=center, scale=scale)
        plot(x=object$cvfit$cv.stats$mean$cv.support,
             y=boxcut.scaled[,1], type='n',
             xlim=range(0,1),
             ylim=range(boxcut.scaled),
             main="Covariate Importance (average values)", cex.main=cex,
             xlab="",
             ylab="", ...)
        for (j in 1:p) {
            lines(x=object$cvfit$cv.stats$mean$cv.support,
                  y=boxcut.scaled[,j],
                  type='l', col=col.cov[j], lty=lty.cov[j], lwd=lwd.cov[j], ...)
        }
        legend("topleft", inset=0.01, legend=varnames[toplot], col=col.cov, lty=lty.cov, lwd=lwd.cov, cex=cex)
        if (center)
            abline(h=0, lty=2, col=1, lwd=0.3, xpd=FALSE)
        if (add.legend)
            legend("bottom", inset=0.01, legend=text.legend, cex=cex)
        mtext(text=xlab, cex=cex, side=1, line=1, outer=FALSE)
        mtext(text=ylab, cex=cex, side=2, line=2, outer=FALSE)

        ticknames <- paste(varnames[toplot], " -", sep="")
        pointtrace <- c(object$cvfit$cv.trace[2], object$cvfit$cv.trace[-1])
        matchtrace <- pmatch(x=pointtrace, table=toplot, duplicates.ok = TRUE)
        plot(x=object$cvfit$cv.stats$mean$cv.support,
             y=matchtrace,
             type='S', yaxt="n", col=col, lty=lty, lwd=lwd,
             xlim=range(0, 1),
             ylim=range(0, p),
             main="Covariate Usage (modal values)", cex.main=cex,
             xlab="",
             ylab="", ...)
        par(mgp=c(1.5, 0, 0))
        axis(side=2, at=1:p, labels=ticknames, tick=FALSE, las=1, line=NA, cex.axis=cex, outer=FALSE)
        if (add.legend)
            legend("bottom", inset=0.01, legend=text.legend, cex=cex)
        mtext(text=xlab, cex=cex, side=1, line=1, outer=FALSE)
        mtext(text="Covariates Used", cex=cex, side=2, line=1+maxlength/2, outer=FALSE)
        if (!is.null(main)) {
            mtext(text=main, cex=1, side=3, outer=TRUE)
        }
    }

    if (is.null(device)) {
        cat("Device: ",  dev.cur(), "\n")
        boxtraceplot(object=object,
                     main=main, xlab=xlab, ylab=ylab,
                     toplot=toplot,
                     center=center, scale=scale,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        cat("Device: ",  dev.cur(), "\n")
        boxtraceplot(object=object,
                     main=main, xlab=xlab, ylab=ylab,
                     toplot=toplot,
                     center=center, scale=scale,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        cat("Device: ",  dev.cur(), "\n")
        boxtraceplot(object=object,
                     main=main, xlab=xlab, ylab=ylab,
                     toplot=toplot,
                     center=center, scale=scale,
                     col.cov=col.cov, lty.cov=lty.cov, lwd.cov=lwd.cov,
                     col=col, lty=lty, lwd=lwd,
                     cex=cex, add.legend=add.legend, text.legend=text.legend)
        dev.off()
    } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
  } else {
    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")
  }
  invisible()
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    plot_boxkm (object,
#                                main=NULL,
#                                xlab="Time", ylab="Probability",
#                                precision=1e-3, mark=3, col=2, cex=1,
#                                steps=1:object$cvfit$cv.nsteps,
#                                nr=3, nc=4,
#                                device=NULL, file="Survival Plots", path=getwd(),
#                                horizontal=TRUE, width=11.5, height=8.5, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

plot_boxkm <- function(object,
                       main=NULL,
                       xlab="Time", ylab="Probability",
                       precision=1e-3, mark=3, col=2, cex=1,
                       steps=1:object$cvfit$cv.nsteps,
                       nr=3, nc=4,
                       device=NULL, file="Survival Plots", path=getwd(),
                       horizontal=TRUE, width=11.5, height=8.5, ...) {

  if (!inherits(object, 'PRSP'))
        stop("Primary argument much be an object of class 'PRSP' \n")

  if (object$plot) {

    boxkmplot <- function(object,
                          main, xlab, ylab,
                          precision, mark, col, cex,
                          steps, nr, nc, ...) {

        if (!is.null(main)) {
            par(mfrow=c(nr, nc), oma=c(0, 0, 3, 0), mar=c(2.5, 2.5, 1.5, 1.5), mgp=c(1.5, 0.5, 0))
        } else {
            par(mfrow=c(nr, nc), oma=c(0, 0, 0, 0), mar=c(2.5, 2.5, 0.0, 1.5), mgp=c(1.5, 0.5, 0))
        }

        times <- object$times
        status <- object$status
        L <- object$cvfit$cv.nsteps
        for (l in steps) {
            boxind <- object$cvfit$cv.boxind[l,]
            ng <- length(unique(boxind[!is.na(boxind)]))
            if (ng == 1) {
                boxind <- 1*boxind
            } else {
                boxind <- 2 - 1*boxind
            }
            surv <- survfit(Surv(times, status) ~ 1 + boxind)
            if (l == 1) {
                plot(surv, main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=2, lwd=0.5, col=col, cex=cex, xlab=xlab, ylab=ylab, ...)
                par(new=TRUE)
                plot(surv, main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=1, lwd=1, col=col, cex=cex, xlab=xlab, ylab=ylab, ...)
            } else {
                plot(surv, main="", conf.int=TRUE, mark.time=FALSE, mark=NA, lty=c(2,2), lwd=c(0.5,0.5), col=c(col,1), cex=cex, xlab=xlab, ylab=ylab, ...)
                par(new=TRUE)
                plot(surv, main="", conf.int=FALSE, mark.time=TRUE, mark=mark, lty=c(1,1), lwd=c(1,1), col=c(col,1), cex=cex, xlab=xlab, ylab=ylab, ...)
            }
            legend("topright", inset=0.01, legend=c("outbox", "inbox"), lty=c(1,1), lwd=c(1,1), col=c(1,2), cex=0.9*cex)
            if (object$cpv) {
                if (object$cvfit$cv.pval[l] <= precision) {
                    legend("bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                           legend=bquote(italic(p) <= .(precision)))
                } else {
                    legend("bottom", inset=0.11, col="black", cex=0.9*cex, bty="n",
                           legend=bquote(italic(p) == .(format(x=object$cvfit$cv.pval[l], scientific=FALSE, digits=4, nsmall=4))))
                }
            }
            legend("bottom", inset=0.01, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LHR) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$cv.lhr[l], digits=3, nsmall=3))))
            legend("bottom", inset=0.06, col="black", cex=0.9*cex, bty="n",
                   legend=substitute(group("", list(paste(italic(LRT) == x, sep="")), ""), list(x=format(x=object$cvfit$cv.stats$mean$cv.lrt[l], digits=3, nsmall=3))))
            legend("bottom", inset=0.16, legend=paste("Step ", l-1, sep=""), col=1, cex=0.9*cex, bty="n")
        }

        if (!is.null(main)) {
            mtext(text=main, cex=cex, side=3, outer=TRUE)
        }
    }

    if (is.null(device)) {
        cat("Device: ",  dev.cur(), "\n")
        boxkmplot(object=object,
                  main=main, xlab=xlab, ylab=ylab,
                  precision=precision, mark=mark, col=col, cex=cex,
                  steps=steps,
                  nr=nr, nc=nc)
    } else if (device == "PS") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".ps", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        postscript(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, horizontal=horizontal)
        cat("Device: ",  dev.cur(), "\n")
        boxkmplot(object=object,
                  main=main, xlab=xlab, ylab=ylab,
                  precision=precision, mark=mark, col=col, cex=cex,
                  steps=steps,
                  nr=nr, nc=nc)
        dev.off()
    } else if (device == "PDF") {
        path <- normalizePath(path=paste(path, "/", sep=""), winslash="\\", mustWork=FALSE)
        file <- paste(file, ".pdf", sep="")
        cat("\nOUTPUT: \n")
        cat("Filename : ", file, "\n")
        cat("Directory: ", path, "\n")
        pdf(file=paste(path, file, sep=""), width=width, height=height, onefile=TRUE, paper=ifelse(test=horizontal, yes="USr", no="US"))
        cat("Device: ",  dev.cur(), "\n")
        boxkmplot(object=object,
                  main=main, xlab=xlab, ylab=ylab,
                  precision=precision, mark=mark, col=col, cex=cex,
                  steps=steps,
                  nr=nr, nc=nc)
        dev.off()
    } else {
        stop("Currently allowed display devices are \"PS\" (Postscript) or \"PDF\" (Portable Document Format) \n")
    }
  } else {
    cat("Either the covariate pre-selection or the Survival Bump Hunting modeling failed for this dataset.\n
        So, there is nothing to plot here.\n")
  }
  invisible()
}
##########################################################################################################################################


