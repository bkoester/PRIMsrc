##########################################################################################################################################
# PRIMsrc
##########################################################################################################################################

##########################################################################################################################################
# SURVIVAL INTERNAL SUBROUTINES
# (never to be called by end-user)
##########################################################################################################################################

##########################################################################################################################################
#################
# Usage         :
################
#                   cv.box.rep(x, times, status,
#                              B, K, arg,
#                              cvtype, 
#                              varsel, varsign, 
#                              probval, timeval,
#                              parallel, seed)
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

cv.box.rep <- function(x, times, status,
                       B, K, arg,
                       cvtype, 
                       varsel, varsign, 
                       probval, timeval,
                       parallel, seed) {

  CV.trace <- vector(mode="list", length=B)
  CV.boxind <- vector(mode="list", length=B)
  CV.boxcut <- vector(mode="list", length=B)
  CV.support <- vector(mode="list", length=B)
  CV.lhr <- vector(mode="list", length=B)
  CV.lrt <- vector(mode="list", length=B)
  CV.cer <- vector(mode="list", length=B)
  CV.time.bar <- vector(mode="list", length=B)
  CV.prob.bar <- vector(mode="list", length=B)
  CV.max.time.bar <- vector(mode="list", length=B)
  CV.min.prob.bar <- vector(mode="list", length=B)
  b <- 1
  k <- 0
  success <- TRUE
  while (b <= B) {
    cat("replicate : ", b, "\n", sep="")
    if (!parallel) {
      set.seed(seed[b])
      cat("seed : ", seed[b], "\n", sep="")
    }
    if (cvtype == "averaged") {
      CVBOX <- cv.ave.box(x=x, times=times, status=status,
                          K=K, arg=arg,
                          varsel=varsel, varsign=varsign, 
                          probval=probval, timeval=timeval,
                          seed=seed[b])
    } else if ((cvtype == "combined") || (cvtype == "none")) {
      CVBOX <- cv.comb.box(x=x, times=times, status=status,
                           K=K, arg=arg,
                           varsel=varsel, varsign=varsign,
                           probval=probval, timeval=timeval,
                           seed=seed[b])
    } else {
      stop("Invalid CV type option \n")
    }
    if (!CVBOX$drop) {
      CV.trace[[b]] <- CVBOX$cvfit$cv.trace
      CV.boxind[[b]] <- CVBOX$cvfit$cv.boxind
      CV.boxcut[[b]] <- CVBOX$cvfit$cv.boxcut
      CV.support[[b]] <- CVBOX$cvfit$cv.stats$cv.support
      CV.lhr[[b]] <- CVBOX$cvfit$cv.stats$cv.lhr
      CV.lrt[[b]] <- CVBOX$cvfit$cv.stats$cv.lrt
      CV.cer[[b]] <- CVBOX$cvfit$cv.stats$cv.cer
      CV.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.time.bar
      CV.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.prob.bar
      CV.max.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.max.time.bar
      CV.min.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.min.prob.bar
      b <- b + 1
      k <- 0
    } else {
      cat("Could not complete one of the cross-validation folds within replicate #", b ,". Retrying replicate with new seed\n", sep="")
      seed[b] <- seed[b] + B
      k <- k + 1
      if (k == 10) {
        cat("Could not complete requested replications after 10 successive trials. Exiting replications\n", sep="")
        b <- B
        success <- FALSE
      }
    }
  }

  return(list("cv.trace"=CV.trace,
              "cv.boxind"=CV.boxind,
              "cv.boxcut"=CV.boxcut,
              "cv.support"=CV.support,
              "cv.lhr"=CV.lhr,
              "cv.lrt"=CV.lrt,
              "cv.cer"=CV.cer,
              "cv.time.bar"=CV.time.bar,
              "cv.prob.bar"=CV.prob.bar,
              "cv.max.time.bar"=CV.max.time.bar,
              "cv.min.prob.bar"=CV.min.prob.bar,
              "success"=success,
              "seed"=seed))
}
##########################################################################################################################################





##########################################################################################################################################
#################
# Usage         :
################
#                   cv.tune.rep(x, times, status,
#                               B, K, arg,
#                               cvtype,
#                               cvcriterion,
#                               fdr, thr,
#                               parallel, seed)
#
################
# Description   :
################
#                    Return internal K-fold (replicated) averaged or combined cross-validated values of
#                    tuned parameter peeling length of box peeling sequence using either the optimal
#                    LHR, LRT or CER as a CV criterion.
#                    CV criterion determined either by:
#                      - maximization of the LHR (between in and out box test samples)
#                      - maximization of the LRT (between in and out box test samples)
#                      - minimization of the CER (between predicted and observed inbox test samples survival times)
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

cv.tune.rep <- function(x, times, status,
                        B, K, arg,
                        cvtype, 
                        cvcriterion,
                        fdr, thr,
                        parallel, seed) {

  # Parsing and evaluating parameters
  alpha <- NULL
  beta <- NULL
  minn <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  # Handling all variable selection cases
  if (is.null(fdr) && is.null(thr)) {
    thr <- ncol(x)
    M <- 1
  } else {
    if (is.null(thr)) {
        M <- length(fdr)
    } else if (is.null(fdr)) {
        M <- length(thr)
    }
  }

  CV.maxsteps <- vector(mode="list", length=B)
  CV.profile <- vector(mode="list", length=B)
  CV.sel <- vector(mode="list", length=B)
  CV.sign <- vector(mode="list", length=B)
  
  b <- 1
  success <- TRUE
  while (b <= B) {

    cat("replicate : ", b, "\n", sep="")
    if (!parallel) {
      set.seed(seed[b])
      cat("seed : ", seed[b], "\n", sep="")
    }
    
    CV.maxsteps.tmp <- numeric(M)
    CV.profile.tmp <- vector(mode="list", length=M)
    CV.sel.tmp <- vector(mode="list", length=M)
    CV.sign.tmp <- vector(mode="list", length=M)

    for (m in 1:M) {
        cat("model : ", m, "\n", sep="")
        if (cvtype == "averaged") {
            CVBOX <- cv.ave.tune(x=x, times=times, status=status,
                                 K=K, alpha=alpha, beta=beta, minn=minn,
                                 fdr=fdr[m], thr=thr[m],
                                 peelcriterion=peelcriterion,
                                 cvcriterion=cvcriterion,
                                 seed=seed[b])
        } else if (cvtype == "combined") {
            CVBOX <- cv.comb.tune(x=x, times=times, status=status,
                                  K=K, alpha=alpha, beta=beta, minn=minn,
                                  fdr=fdr[m], thr=thr[m],
                                  peelcriterion=peelcriterion,
                                  cvcriterion=cvcriterion,
                                  seed=seed[b])
        } else if (cvtype == "none") {
            CVBOX <- cv.comb.tune(x=x, times=times, status=status,
                                  K=1, alpha=alpha, beta=beta, minn=minn,
                                  fdr=fdr, thr=thr,
                                  peelcriterion=peelcriterion,
                                  cvcriterion=NULL,
                                  seed=seed[b])
        } else {
            stop("Invalid CV type option \n")
        }
        if (!CVBOX$drop) {
            CV.maxsteps.tmp[m] <- CVBOX$cvfit$cv.maxsteps
            CV.profile.tmp[[m]] <- CVBOX$cvfit$cv.profile
            CV.sel.tmp[[m]] <- CVBOX$cvfit$cv.sel
            CV.sign.tmp[[m]] <- CVBOX$cvfit$cv.sign
        } else {
            k <- 0
            retry <- TRUE
            while (retry) {
                cat("Warning: could not complete one of the cross-validation folds for model #", m, " within replicate #", b, ". 
                    Retrying replicate for that model with new seed\n", sep="")
                seed[b] <- seed[b] + B
                if (cvtype == "averaged") {
                    CVBOX <- cv.ave.tune(x=x, times=times, status=status,
                                         K=K, alpha=alpha, beta=beta, minn=minn,
                                         fdr=fdr[m], thr=thr[m],
                                         peelcriterion=peelcriterion,
                                         cvcriterion=cvcriterion,
                                         seed=seed[b])
                } else if (cvtype == "combined") {
                    CVBOX <- cv.comb.tune(x=x, times=times, status=status,
                                          K=K, alpha=alpha, beta=beta, minn=minn,
                                          fdr=fdr[m], thr=thr[m],
                                          peelcriterion=peelcriterion,
                                          cvcriterion=cvcriterion,
                                          seed=seed[b])
                } else if (cvtype == "none") {
                    CVBOX <- cv.comb.tune(x=x, times=times, status=status,
                                          K=1, alpha=alpha, beta=beta, minn=minn,
                                          fdr=fdr, thr=thr,
                                          peelcriterion=peelcriterion,
                                          cvcriterion=NULL,
                                          seed=seed[b])
                } else {
                    stop("Invalid CV type option \n")
                }
                if (!CVBOX$drop) {
                    CV.maxsteps.tmp[m] <- CVBOX$cvfit$cv.maxsteps
                    CV.profile.tmp[[m]] <- CVBOX$cvfit$cv.profile
                    CV.sel.tmp[[m]] <- CVBOX$cvfit$cv.sel
                    CV.sign.tmp[[m]] <- CVBOX$cvfit$cv.sign
                    retry <- FALSE
                } else {
                    k <- k + 1            
                }
            }
            if (k == 10) {
                cat("Could not complete requested replications after 10 successive trials. Exiting replications\n", sep="")
                b <- B
                success <- FALSE
            }
        }
    }
    CV.maxsteps[[b]] <- ceiling(mean(CV.maxsteps.tmp))
    if (cvtype == "none") {
        CV.profile[[b]] <- NULL
    } else if ((cvtype == "averaged") || (cvtype == "combined")) {
        if ((cvcriterion == "lhr") || (cvcriterion == "lrt")) {
                CV.profile[[b]] <- list2mat(list=CV.profile.tmp, fill=0, coltrunc=CV.maxsteps[[b]])
        } else if (cvcriterion == "cer") {
                CV.profile[[b]] <- list2mat(list=CV.profile.tmp, fill=1, coltrunc=CV.maxsteps[[b]])
        }
        dimnames(CV.profile[[b]]) <- list(paste("model", 1:M, sep=""), paste("step", 0:(CV.maxsteps[[b]]-1), sep=""))
    } else {
        stop("Invalid CV type option \n")
    }
    CV.sel[[b]] <- CV.sel.tmp
    CV.sign[[b]] <- CV.sign.tmp
    b <- b + 1

  }

  return(list("cv.maxsteps"=CV.maxsteps,
              "cv.profile"=CV.profile,
              "cv.sel"=CV.sel,
              "cv.sign"=CV.sign,
              "success"=success,
              "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.pval (x, times, status,
#                             cvtype,
#                             varsel, varsign,
#                             A, K, arg, 
#                             obs.chisq,
#                             parallel, conf)
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

cv.pval <- function(x, times, status,
                    cvtype,
                    varsel, varsign,
                    A, K, arg, 
                    obs.chisq,
                    parallel, conf) {

  if (!parallel) {
    null.chisq <- cv.null(x=x, times=times, status=status,
                          cvtype=cvtype,
                          varsel=varsel, varsign=varsign,
                          A=A, K=K, arg=arg)
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
    clusterSetRNGStream(cl=cl, iseed=NULL)
    null.cl <- clusterCall(cl=cl, fun=cv.null,
                           x=x, times=times, status=status,
                           cvtype=cvtype,
                           varsel=varsel, varsign=varsign,
                           A=ceiling(A/conf$cpus), K=K, arg=arg)
    stopCluster(cl)
    null.chisq <- cbindlist(null.cl)
  }
  cvl <- nrow(null.chisq)
  pval <- numeric(cvl)
  for (l in 1:cvl) {
    pval[l] <- mean((null.chisq[l,] >= obs.chisq[l]), na.rm=TRUE)
  }

  return(pval)
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.null (x, times, status,
#                            cvtype,
#                            varsel, varsign,
#                            A, K, arg)
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

cv.null <- function(x, times, status,
                    cvtype,
                    varsel, varsign,
                    A, K, arg) {

  n <- nrow(x)
  null.chisq <- vector(mode="list", length=A)
  a <- 1
  while (a <= A) {
    cat("Permutation sample: ", a, "\n")
    perm.ind <- sample(x = 1:n, size = n, replace = FALSE, prob = NULL)
    perm.times <- times[perm.ind]
    perm.status <- status[perm.ind]
    if (cvtype == "averaged") {
      obj <- tryCatch({cv.ave.box(x=x, times=perm.times, status=perm.status, 
                                  varsel=varsel, varsign=varsign,
                                  probval=NULL, timeval=NULL, 
                                  K=K, arg=arg, seed=NULL)
                      }, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else if ((cvtype == "combined") || (cvtype == "none")) {
      obj <- tryCatch({cv.comb.box(x=x, times=perm.times, status=perm.status, 
                                   varsel=varsel, varsign=varsign,
                                   probval=NULL, timeval=NULL, 
                                   K=K, arg=arg, seed=NULL)
                      }, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else {
      stop("Invalid CV type option \n")
    }
  }

  return(t(list2mat(null.chisq)))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.ave.tune(x, times, status, 
#                               K, alpha, beta, minn,
#                               fdr, thr,
#                               peelcriterion,
#                               cvcriterion,
#                               seed)
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
                            
cv.ave.tune <- function(x, times, status, 
                        K, alpha, beta, minn, 
                        fdr, thr,
                        peelcriterion, 
                        cvcriterion,
                        seed) {

  folds <- cv.folds(n=nrow(x), K=K, seed=seed)
  maxsteps <- numeric(K)
  sel.list <- vector(mode="list", length=K)
  sign.list <- vector(mode="list", length=K)
  boxstat.list <- vector(mode="list", length=K)

  k <- 1
  while (k <= K) {
    cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }

    # Training the model
    peelobj <- cv.tune(traindata=traindata, traintime=traintime, trainstatus=trainstatus, 
                       alpha=alpha, beta=beta, minn=minn, 
                       fdr=fdr, thr=thr,
                       peelcriterion=peelcriterion,
                       seed=seed)
    if (is.empty(peelobj$varsel)) {
        drop <- TRUE
        k <- K + 1
    } else {
        maxsteps[k] <- peelobj$maxsteps
        testdata <- testdata[,peelobj$varsel, drop=FALSE]
        # Compute the box statistics for all steps, each entry or row signifies a step
        boxstat <- vector(mode="list", length=maxsteps[k])
        ind.rem <- numeric(0)
        for (l in 1:maxsteps[k]) {
            # Extract the rule and sign as one vector
            boxcut <- peelobj$boxcut[l,] * peelobj$varsign
            test.cut <- t(t(testdata) * peelobj$varsign)
            test.ind <- t(t(test.cut) >= boxcut)
            # Set as TRUE which observations are TRUE for all covariates
            test.ind <- (rowMeans(test.ind) == 1)
            test.ind1 <- 1*test.ind
            if ((l == 1) && (sum(test.ind, na.rm=TRUE) != 0)) {
                if (is.null(cvcriterion)) {
                    boxstat[[l]] <- NA
                } else if (cvcriterion == "lhr") {
                    boxstat[[l]] <- 0
                } else if (cvcriterion == "lrt") {
                    boxstat[[l]] <- 0
                } else if (cvcriterion == "cer") {
                    boxstat[[l]] <- 1
                }
            } else if ((sum(test.ind, na.rm=TRUE) != length(test.ind[!is.na(test.ind)])) && (sum(test.ind, na.rm=TRUE) != 0)) {
                surv.formula <- (Surv(testtime, teststatus) ~ 1 + test.ind1)
                if (is.null(cvcriterion)) {
                    boxstat[[l]] <- NA
                } else if (cvcriterion == "lhr") {
                    coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
                    boxstat[[l]] <- coxobj$coef
                } else if (cvcriterion == "lrt") {
                    boxstat[[l]] <- survdiff(surv.formula, rho=0)$chisq
                } else if (cvcriterion == "cer") {
                    coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
                    predobj <- predict(object=coxobj, type="lp", reference="sample")
                    boxstat[[l]] <- rcorr.cens(x=predobj, S=Surv(testtime, teststatus))['C Index']
                }
            } else {
                if (is.null(cvcriterion)) {
                    boxstat[[l]] <- NA
                } else if (cvcriterion == "lhr") {
                    boxstat[[l]] <- 0
                } else if (cvcriterion == "lrt") {
                    boxstat[[l]] <- 0
                } else if (cvcriterion == "cer") {
                    boxstat[[l]] <- 1
                }
                ind.rem <- c(ind.rem, l)
            }
        }
        # Store the test set results from each fold
        boxstat.list[[k]] <- boxstat
        sel.list[[k]] <- peelobj$varsel
        sign.list[[k]] <- peelobj$varsign
        # Drop case
        if (length(ind.rem) == maxsteps[k]) {
            drop <- TRUE
            k <- K + 1
        } else {
            drop <- FALSE
            k <- k + 1
        }
    }
  }

  if (drop) {     
  
    CV.Lm <- 0
    CV.profile <- NULL
    CV.sel <- NULL
    CV.sign <- NULL

  } else {
    
    # Cross-validated minimum length from all folds
    CV.Lm <- min(maxsteps)
    
    # Get the selected covariates from all folds
    varfreq <- table(unlist(sel.list), useNA="no")
    #w <- varfreq[varfreq == K]           # most conservative (any covariate present in all the folds)
    w <- varfreq[varfreq >= ceiling(K/2)] # medium conservative (any covariate present in at least half of the folds)
    #w <- varfreq                         # least conservative (any covariate present in at least one of the folds)
    CV.sel <- as.numeric(names(w))
    names(CV.sel) <- colnames(x)[CV.sel]
    
    # Get the covariates directions of peeling from all folds
    CV.sign <- sapply(X=sign.list, FUN=function(x, id){x[id]}, id=names(CV.sel), simplify=TRUE)
    if (is.vector(CV.sign)) CV.sign <- matrix(data=CV.sign, nrow=1, ncol=K, dimnames=list(names(CV.sel),NULL))
    CV.sign <- apply(X=CV.sign, MARGIN=1, FUN=function(x){as.numeric(names(which.max(table(x, useNA="no"))))})
    names(CV.sign) <- names(CV.sel)

    # Truncate the cross-validated quantities from all folds to the same cross-validated length
    for (k in 1:K) {
        boxstat.list[[k]] <- boxstat.list[[k]][1:CV.Lm]
    }
    
    # Compute the averaged box statistics for each step from all folds
    # Each entry or row signifies a step
    CV.profile <- numeric(CV.Lm)
    names(CV.profile) <- paste("step", 0:(CV.Lm-1), sep="")
    if (is.null(cvcriterion)) {
        CV.profile <- rep(NA, CV.Lm)
    } else if (cvcriterion == "lhr") {
        for (l in 1:CV.Lm) {
            sumlhr <- rep(NA, K)
            for (k in 1:K) {
                sumlhr[k] <- boxstat.list[[k]][[l]]
            }
            CV.profile[l] <- mean(sumlhr, na.rm=TRUE)
        }
    } else if (cvcriterion == "lrt") {
        for (l in 1:CV.Lm) {
            sumlrt <- rep(NA, K)
            for (k in 1:K) {
                sumlrt[k] <- boxstat.list[[k]][[l]]
            }
            CV.profile[l] <- mean(sumlrt, na.rm=TRUE)
        }
    } else if (cvcriterion == "cer") {
        for (l in 1:CV.Lm) {
            sumcer <- rep(NA, K)
            for (k in 1:K) {
                sumcer[k] <- boxstat.list[[k]][[l]]
            }
            CV.profile[l] <- mean(sumcer, na.rm=TRUE)
        }
    }
    
  }

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.Lm,
                 "cv.profile"=CV.profile,
                 "cv.sel"=CV.sel,
                 "cv.sign"=CV.sign)

  return(list("cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################





##########################################################################################################################################
#################
# Usage         :
################
#                   cv.comb.tune(x, times, status, 
#                                K, alpha, beta, minn,
#                                fdr, thr,
#                                peelcriterion,
#                                cvcriterion,
#                                seed)
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

cv.comb.tune <- function(x, times, status, 
                         K, alpha, beta, minn, 
                         fdr, thr,
                         peelcriterion, 
                         cvcriterion,
                         seed) {

  folds <- cv.folds(n=nrow(x), K=K, seed=seed)
  maxsteps <- numeric(K)
  times.list <- vector(mode="list", length=K)
  status.list <- vector(mode="list", length=K)
  sel.list <- vector(mode="list", length=K)
  sign.list <- vector(mode="list", length=K)
  boxind.list <- vector(mode="list", length=K)

  k <- 1
  while (k <= K) {
    cat("Fold : ", k, "\n", sep="")
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }

    # Store the test set data from each fold (Note: the order of observations of times and status from each fold is kept in the list)
    times.list[[k]] <- testtime
    status.list[[k]] <- teststatus
    
    # Training the model
    peelobj <- cv.tune(traindata=traindata, traintime=traintime, trainstatus=trainstatus, 
                       alpha=alpha, beta=beta, minn=minn, 
                       fdr=fdr, thr=thr,
                       peelcriterion=peelcriterion,
                       seed=seed)
    if (is.empty(peelobj$varsel)) {
        drop <- TRUE
        k <- K + 1
    } else {
        maxsteps[k] <- peelobj$maxsteps
        testdata <- testdata[, peelobj$varsel, drop=FALSE]
        # Create the indicator matrix of the test data that is within the box for each step
        boxind <- matrix(NA, nrow=maxsteps[k], ncol=nrow(testdata))
        for (l in 1:maxsteps[k]) {
            # Extract the rule and sign as one vector
            boxcut <- peelobj$boxcut[l,] * peelobj$varsign
            test.cut <- t(t(testdata) * peelobj$varsign)
            test.ind <- t(t(test.cut) >= boxcut)
            # Set as TRUE which observations are TRUE for all covariates
            boxind[l, ] <- (rowMeans(test.ind) == 1)
        }
        # Store the test set results from each fold
        sel.list[[k]] <- peelobj$varsel
        sign.list[[k]] <- peelobj$varsign
        boxind.list[[k]] <- boxind
        drop <- FALSE
        k <- k + 1
    }
  }

  if (drop) {
    
    CV.Lm <- 0
    CV.profile <- NULL
    CV.sel <- NULL
    CV.sign <- NULL
    
  } else {
  
    # Cross-validated minimum length from all folds
    CV.Lm <- min(maxsteps)
    
    # Concatenates the observations of test times and status from all folds
    # Re-ordered by initial order of observations
    ord <- folds$key
    CV.times <- unlist(times.list)[ord]
    CV.status <- unlist(status.list)[ord]
    
    # Get the selected covariates from all folds
    varfreq <- table(unlist(sel.list), useNA="no")
    #w <- varfreq[varfreq == K]           # most conservative (any covariate present in all the folds)
    w <- varfreq[varfreq >= ceiling(K/2)] # medium conservative (any covariate present in at least half of the folds)
    #w <- varfreq                         # least conservative (any covariate present in at least one of the folds)
    CV.sel <- as.numeric(names(w))
    names(CV.sel) <- colnames(x)[CV.sel]
    
    # Get the covariates directions of peeling from all folds
    CV.sign <- sapply(X=sign.list, FUN=function(x, id){x[id]}, id=names(CV.sel), simplify=TRUE)
    if (is.vector(CV.sign)) CV.sign <- matrix(data=CV.sign, nrow=1, ncol=K, dimnames=list(names(CV.sel),NULL))
    CV.sign <- apply(X=CV.sign, MARGIN=1, FUN=function(x){as.numeric(names(which.max(table(x, useNA="no"))))})
    names(CV.sign) <- names(CV.sel)
    
    # Get the test box membership indicator vector of all observations for each step from all folds
    # Based on the combined membership indicator vectors over the folds
    # Re-ordered by initial order of observations
    CV.boxind <- cbindlist(boxind.list, trunc=CV.Lm)[,ord]
    rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
    colnames(CV.boxind) <- rownames(x)
     
    # General case: compute the combined test box statistics from all folds for all steps, each entry or row signifies a step
    CV.profile <- numeric(CV.Lm)
    names(CV.profile) <- paste("step", 0:(CV.Lm-1), sep="")
    ind.rem <- numeric(0)
    if (is.null(cvcriterion)) {
        CV.profile <- rep(NA, CV.Lm)
    } else if (cvcriterion == "lhr") {
        for (l in 1:CV.Lm) {
            boxind <- CV.boxind[l,]
            boxind1 <- 1*boxind
            if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
                CV.profile[l] <- 0
            } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
                surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
                coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
                CV.profile[l] <- coxobj$coef
            } else {
                CV.profile[l] <- 0
                ind.rem <- c(ind.rem, l)
            }
        }
    } else if (cvcriterion == "lrt") {
        for (l in 1:CV.Lm) {
            boxind <- CV.boxind[l,]
            boxind1 <- 1*boxind
            if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
                CV.profile[l] <- 0
            } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
                surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
                CV.profile[l] <- survdiff(surv.formula, rho=0)$chisq
            } else {
                CV.profile[l] <- 0
                ind.rem <- c(ind.rem, l)
            }
        }
    } else if (cvcriterion == "cer") {
        for (l in 1:CV.Lm) {
            boxind <- CV.boxind[l,]
            boxind1 <- 1*boxind
            if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
                CV.profile[l] <- 1
            } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
                surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
                coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
                predobj <- predict(object=coxobj, type="lp", reference="sample")
                CV.profile[l] <- rcorr.cens(x=predobj, S=Surv(CV.times, CV.status))['C Index']
            } else {
                CV.profile[l] <- 1
                ind.rem <- c(ind.rem, l)
            }
        }
    }
    
    # Drop case
    if (length(ind.rem) == CV.Lm) {
        drop <- TRUE
        CV.Lm <- NULL
        CV.profile <- NULL
        CV.sel <- NULL
        CV.sign <- NULL
    }
    
  }

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.Lm,
                 "cv.profile"=CV.profile,
                 "cv.sel"=CV.sel,
                 "cv.sign"=CV.sign)

  return(list("cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.ave.box (x, times, status,
#                               varsel, varsign,
#                               probval, timeval,
#                               K, arg, seed)
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

cv.ave.box <- function(x, times, status,
                       varsel, varsign,
                       probval, timeval,
                       K, arg, seed) {

  fold.obj <- cv.ave.fold(x=x, times=times, status=status,
                          varsign=varsign,
                          probval=probval, timeval=timeval,
                          K=K, arg=arg, seed=seed)
  n <- nrow(x)
  p <- ncol(x)
  
  drop <- fold.obj$drop
  if (drop) {
  
    CV.Lm <- 0
    CV.trace <- NULL
    CV.boxind <- NULL
    CV.boxcut <- NULL
    CV.rules <- NULL
    CV.stats <- NULL
    
  } else {
  
    trace.list <- fold.obj$trace.list
    boxcut.list <- fold.obj$boxcut.list
    boxstat.list <- fold.obj$boxstat.list

    # Cross-validated minimum length from all folds
    CV.Lm <- min(fold.obj$maxsteps)

    # Get the covariates traces from all folds
    # Covariates traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
    CV.trace <- list2mat(list=trace.list, coltrunc=CV.Lm)
    CV.trace <- t(CV.trace)
    dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

    # Truncate the cross-validated quantities from all folds to the same cross-validated length
    for (k in 1:K) {
        boxcut.list[[k]] <- boxcut.list[[k]][1:CV.Lm,names(varsel)]
        boxstat.list[[k]] <- boxstat.list[[k]][1:CV.Lm]
    }

    # Compute the averaged box statistics for each step from all folds
    # Each entry or row signifies a step
    CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), names(varsel)))
    CV.support <- rep(NA, CV.Lm)
    names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.lhr <- rep(NA, CV.Lm)
    names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.lrt <- rep(NA, CV.Lm)
    names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.cer <- rep(NA, CV.Lm)
    names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.time.bar <- rep(NA, CV.Lm)
    names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.prob.bar <- rep(NA, CV.Lm)
    names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.max.time.bar <- rep(NA, CV.Lm)
    names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
    CV.min.prob.bar <- rep(NA, CV.Lm)
    names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
    for (l in 1:CV.Lm) {
        summincut <- matrix(NA, K, p)
        sumtime <- rep(NA, K)
        sumprob <- rep(NA, K)
        summaxtime <- rep(NA, K)
        summinprob <- rep(NA, K)
        sumlhr <- rep(NA, K)
        sumlrt <- rep(NA, K)
        sumcer <- rep(NA, K)
        sumsupport <- rep(NA, K)
        for (k in 1:K) {
            summincut[k,] <- boxcut.list[[k]][l,]
            sumlhr[k] <- boxstat.list[[k]][[l]][[1]]
            sumlrt[k] <- boxstat.list[[k]][[l]][[2]]
            sumcer[k] <- boxstat.list[[k]][[l]][[3]]
            sumsupport[k] <- boxstat.list[[k]][[l]][[4]]
            sumtime[k] <- boxstat.list[[k]][[l]][[5]]
            sumprob[k] <- boxstat.list[[k]][[l]][[6]]
            summaxtime[k] <- boxstat.list[[k]][[l]][[7]]
            summinprob[k] <- boxstat.list[[k]][[l]][[8]]
        }
        CV.boxcut[l, ] <- colMeans(summincut, na.rm=TRUE)
        CV.lhr[l] <- mean(sumlhr, na.rm=TRUE)
        CV.lrt[l] <- mean(sumlrt, na.rm=TRUE)
        CV.cer[l] <- mean(sumcer, na.rm=TRUE)
        CV.support[l] <- mean(sumsupport, na.rm=TRUE)
        CV.time.bar[l] <- mean(sumtime, na.rm=TRUE)
        CV.prob.bar[l] <- mean(sumprob, na.rm=TRUE)
        CV.max.time.bar[l] <- mean(summaxtime, na.rm=TRUE)
        CV.min.prob.bar[l] <- mean(summinprob, na.rm=TRUE)
    }

    # Get the box membership indicator vector of all observations for each step from all folds
    # Based on the corresponding averaged box over the folds
    CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
    for (l in 1:CV.Lm) {
        boxcut <- CV.boxcut[l,] * varsign
        x.cut <- t(t(x) * varsign)
        x.ind <- t(t(x.cut) >= boxcut)
        CV.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boudaries for all axes directions
    }
    rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
    colnames(CV.boxind) <- rownames(x)
    
    # Box peeling rules for each step from all folds
    CV.rules <- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), names(varsel))))
    for (l in 1:CV.Lm) {
        for (j in 1:p) {
            if (varsign[j] > 0) {
                ss <- ">="
            } else {
                ss <- "<="
            }
            CV.rules[l,j] <- paste(names(varsel)[j], ss, format(x=CV.boxcut[l,j], digits=3, nsmall=3), sep="")
        }
    }

    # Box statistics for each step
    CV.stats <-  data.frame("cv.support"=CV.support,
                            "cv.lhr"=CV.lhr,
                            "cv.lrt"=CV.lrt,
                            "cv.cer"=CV.cer,
                            "cv.time.bar"=CV.time.bar,
                            "cv.prob.bar"=CV.prob.bar,
                            "cv.max.time.bar"=CV.max.time.bar,
                            "cv.min.prob.bar"=CV.min.prob.bar)
    rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")
  }

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.Lm,
                 "cv.boxind"=CV.boxind,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                   cv.comb.box (x, times, status,
#                                varsel, varsign,
#                                probval, timeval,
#                                K, arg, seed)
#
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

cv.comb.box <- function(x, times, status,
                        varsel, varsign,
                        probval, timeval,
                        K, arg, seed) {

  fold.obj <- cv.comb.fold(x=x, times=times, status=status,
                           varsign=varsign,
                           K=K, arg=arg, seed=seed)
  n <- nrow(x)
  p <- ncol(x)
  cvtimes.list <- fold.obj$cvtimes.list
  cvstatus.list <- fold.obj$cvstatus.list
  trace.list <- fold.obj$trace.list
  boxind.list <- fold.obj$boxind.list
  boxcut.list <- fold.obj$boxcut.list
                  
  # Cross-validated minimum length from all folds
  CV.Lm <- min(fold.obj$maxsteps)

  # Concatenates the observations of test times and status from all folds
  # Re-ordered by initial order of observations
  ord <- fold.obj$key
  CV.times <- unlist(cvtimes.list)[ord]
  CV.status <- unlist(cvstatus.list)[ord]

  # Get the covariates traces from all folds
  # Covariates traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- list2mat(list=trace.list, coltrunc=CV.Lm)
  CV.trace <- t(CV.trace)
  dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

  # Get the test box membership indicator vector of all observations for each step from all folds
  # Based on the combined membership indicator vectors over the folds
  # Re-ordered by initial order of observations
  CV.boxind <- cbindlist(boxind.list, trunc=CV.Lm)[,ord]
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(x)

  # Get the combined boxcut (truncated to the same cross-validated length) for each step from all folds
  # using the circumscribing box to the conmbined test set in-box samples over all the folds
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), names(varsel)))
  tmparray <- list2array(list=boxcut.list, rowtrunc=CV.Lm)
  for (l in 1:CV.Lm) {
        for (j in 1:p) {
            if (varsign[j] > 0) {
                CV.boxcut[l,j] <- min(x[CV.boxind[l,],j])
            } else {
                CV.boxcut[l,j] <- max(x[CV.boxind[l,],j])
            }
        }
  }

  # Box peeling rules for each step from all folds
  CV.rules <- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), names(varsel))))
  for (l in 1:CV.Lm) {
        for (j in 1:p) {
            if (varsign[j] > 0) {
                ss <- ">="
            } else {
                ss <- "<="
            }
            CV.rules[l,j] <- paste(names(varsel)[j], ss, format(x=CV.boxcut[l,j], digits=3, nsmall=3), sep="")
        }
  }

  # Compute the combined test box statistics from all folds for all steps, each entry or row signifies a step
  CV.support <- rep(NA, CV.Lm)
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lhr <- rep(NA, CV.Lm)
  names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lrt <- rep(NA, CV.Lm)
  names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.cer <- rep(NA, CV.Lm)
  names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(NA, CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(NA, CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(NA, CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(NA, CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  timemat <- matrix(NA, nrow=CV.Lm, ncol=n)
  probmat <- matrix(NA, nrow=CV.Lm, ncol=n)
  ind.rem <- numeric(0)
  for (l in 1:CV.Lm) {
        boxind <- CV.boxind[l,]
        boxind1 <- 1*boxind
        if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
            surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
            timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
            probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
            CV.lhr[l] <- 0
            CV.lrt[l] <- 0
            CV.cer[l] <- 1
            CV.support[l] <- 1
        } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
            surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
            timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
            probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
            surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
            coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
            CV.lhr[l] <- coxobj$coef
            CV.lrt[l] <- survdiff(surv.formula, rho=0)$chisq
            predobj <- predict(object=coxobj, type="lp", reference="sample")
            CV.cer[l] <- rcorr.cens(x=predobj, S=Surv(CV.times, CV.status))['C Index']
            CV.support[l] <- mean(boxind, na.rm=TRUE)
        } else {
            timemat[l, ] <- NA
            probmat[l, ] <- NA
            CV.lhr[l] <- 0
            CV.lrt[l] <- 0
            CV.cer[l] <- 1
            CV.support[l] <- NA
            ind.rem <- c(ind.rem, l)
        }
  }
  
  if (length(ind.rem) != CV.Lm) {
  
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:CV.Lm) {
        if (!(l %in% ind.rem)) {
            CV.time.bar[l] <- time.bar[l]
            CV.prob.bar[l] <- prob.bar[l]
            CV.max.time.bar[l] <-  max.time.bar[l]
            CV.min.prob.bar[l] <- min.prob.bar[l]
        } else {
            CV.time.bar[l] <- NA
            CV.prob.bar[l] <- NA
            CV.max.time.bar[l] <- NA
            CV.min.prob.bar[l] <- NA
        }
    }
    CV.stats <- data.frame("cv.support"=CV.support,
                           "cv.lhr"=CV.lhr,
                           "cv.lrt"=CV.lrt,
                           "cv.cer"=CV.cer,
                           "cv.time.bar"=CV.time.bar,
                           "cv.prob.bar"=CV.prob.bar,
                           "cv.max.time.bar"=CV.max.time.bar,
                           "cv.min.prob.bar"=CV.min.prob.bar)
    rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")
        
  } else {
  
    drop <- TRUE
    CV.Lm <- 0
    CV.trace <- NULL
    CV.boxind <- NULL
    CV.boxcut <- NULL
    CV.rules <- NULL
    CV.stats <- NULL
        
  }
    
  # Create the return object 'CV.fit'
  CV.fit <- list("cv.maxsteps"=CV.Lm,
                 "cv.boxind"=CV.boxind,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.ave.fold (x, times, status,
#                                 varsign,
#                                 probval, timeval,
#                                 K, arg, seed)
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

cv.ave.fold <- function(x, times, status,
                        varsign,
                        probval, timeval,
                        K, arg, seed) {

  drop <- FALSE
  folds <- cv.folds(n=nrow(x), K=K, seed=seed)

  maxsteps <- numeric(K)
  boxstat.list <- vector(mode="list", length=K)
  boxcut.list <- vector(mode="list", length=K)
  trace.list <- vector(mode="list", length=K)

  for (k in 1:K) {
    cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }
    peelobj <- cv.ave.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                           testdata=testdata, teststatus=teststatus, testtime=testtime,
                           varsign=varsign,
                           probval=probval, timeval=timeval,
                           arg=arg, seed=seed)

    # Store the test set results from each fold
    maxsteps[k] <- peelobj$maxsteps
    boxstat.list[[k]] <- peelobj$boxstat
    boxcut.list[[k]] <- peelobj$boxcut
    trace.list[[k]] <- peelobj$vartrace
    drop <- (drop || peelobj$drop)
  }

  return(list("maxsteps"=maxsteps,
              "boxstat.list"=boxstat.list,
              "boxcut.list"=boxcut.list,
              "trace.list"=trace.list,
              "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.comb.fold (x, times, status,
#                                  varsign,
#                                  K, arg, seed)
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

cv.comb.fold <- function(x, times, status,
                         varsign,
                         K, arg, seed) {

  folds <- cv.folds(n=nrow(x), K=K, seed=seed)
  maxsteps <- numeric(K)
  cvtimes.list <- vector(mode="list", length=K)
  cvstatus.list <- vector(mode="list", length=K)
  boxind.list <- vector(mode="list", length=K)
  boxcut.list <- vector(mode="list", length=K)
  trace.list <- vector(mode="list", length=K)

  k <- 1
  while (k <= K) {
    cat("Fold : ", k, "\n", sep="")
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }

    # Store the test set data from each fold (Note: the order of observations of times and status from each fold is kept in the list)
    cvtimes.list[[k]] <- testtime
    cvstatus.list[[k]] <- teststatus

    peelobj <- cv.comb.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                            testdata=testdata, teststatus=teststatus, testtime=testtime,
                            varsign=varsign,
                            arg=arg, seed=seed)

    # Store the test set results from each fold
    maxsteps[k] <- peelobj$maxsteps
    boxind.list[[k]] <- peelobj$boxind
    boxcut.list[[k]] <- peelobj$boxcut
    trace.list[[k]] <- peelobj$vartrace
    k <- k + 1
  }

  return(list("maxsteps"=maxsteps,
              "key"=folds$key,
              "cvtimes.list"=cvtimes.list,
              "cvstatus.list"=cvstatus.list,
              "boxcut.list"=boxcut.list,
              "boxind.list"=boxind.list,
              "trace.list"=trace.list))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.ave.peel (traindata, trainstatus, traintime,
#                                 testdata, teststatus, testtime,
#                                 varsign,
#                                 probval, timeval,
#                                 arg, seed)
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

cv.ave.peel <- function(traindata, trainstatus, traintime,
                        testdata, teststatus, testtime,
                        varsign,
                        probval, timeval,
                        arg, seed) {

  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus, varsign=varsign, arg=arg, seed=seed)
  maxsteps <- peelobj$maxsteps

  # Compute the box statistics for all steps, each entry or row signifies a step
  boxstat <- vector(mode="list", length=maxsteps)
  timemat <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  probmat <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  ind.rem <- numeric(0)
  for (l in 1:maxsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l,] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    test.ind <- (rowMeans(test.ind) == 1)
    test.ind1 <- 1*test.ind
    if ((l == 1) && (sum(test.ind, na.rm=TRUE) != 0)) {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- 1
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else if ((sum(test.ind, na.rm=TRUE) != length(test.ind[!is.na(test.ind)])) && (sum(test.ind, na.rm=TRUE) != 0)) {
      surv.formula <- (Surv(testtime, teststatus) ~ 1 + test.ind1)
      coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
      lhr <- coxobj$coef
      lrt <- survdiff(surv.formula, rho=0)$chisq
      predobj <- predict(object=coxobj, type="lp", reference="sample")
      cer <- rcorr.cens(x=predobj, S=Surv(testtime, teststatus))['C Index']
      support <- mean(test.ind)
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- NA
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      timemat[l, ] <- NA
      probmat[l, ] <- NA
      ind.rem <- c(ind.rem, l)
    }
  }
  if (length(ind.rem) != maxsteps) {
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:maxsteps) {
      if (!(l %in% ind.rem)) {
        boxstat[[l]] <- c(boxstat[[l]], time.bar[l], prob.bar[l], max.time.bar[l], min.prob.bar[l])
      } else {
        boxstat[[l]] <- c(boxstat[[l]], NA, NA, NA, NA)
      }
    }
  } else {
    drop <- TRUE
    for (l in 1:maxsteps) {
      boxstat[[l]] <- rep(NA, 8)
    }
  }
  
  return(list("maxsteps"=maxsteps, 
              "boxcut"=peelobj$boxcut, 
              "boxstat"=boxstat, 
              "vartrace"=peelobj$vartrace, 
              "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.comb.peel (traindata, trainstatus, traintime,
#                                  testdata, teststatus, testtime,
#                                  varsign,
#                                  arg, seed)
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

cv.comb.peel <- function(traindata, trainstatus, traintime,
                         testdata, teststatus, testtime,
                         varsign,
                         arg, seed) {

  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus, varsign=varsign, arg=arg, seed=seed)
  maxsteps <- peelobj$maxsteps
  
  # Create the indicator matrix of the test data that is within the box for each step
  boxind <- matrix(NA, nrow=maxsteps, ncol=nrow(testdata))
  for (l in 1:maxsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l,] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    boxind[l, ] <- (rowMeans(test.ind) == 1)
  }      

  return(list("maxsteps"=maxsteps, 
              "boxcut"=peelobj$boxcut, 
              "boxind"=boxind, 
              "vartrace"=peelobj$vartrace))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    peel.box (traindata, traintime, trainstatus, varsign, arg, seed)
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

peel.box <- function(traindata, traintime, trainstatus, varsign, arg, seed) {

  if (!is.null(seed))
    set.seed(seed)

  # Parsing and evaluating parameters
  alpha <- NULL
  beta <- NULL
  minn <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  # Constants
  n <- nrow(traindata)                                   # Number of samples
  p <- ncol(traindata)                                   # Number of initially pre-selected covariates
  beta <- max(minn/n, beta)                              # Minimal box support

  # Fitting by the PRSP algorithm
  prsp.fit <- prsp(traindata=traindata, 
                   traintime=traintime, 
                   trainstatus=trainstatus,
                   L=L,
                   varsign=varsign, 
                   alpha=alpha, 
                   beta=beta,
                   peelcriterion=peelcriterion)

  # Returning the final results
  return(list("maxsteps"=prsp.fit$maxsteps,
              "vartrace"=prsp.fit$vartrace,
              "boxcut"=prsp.fit$boxcut))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.tune(traindata, traintime, trainstatus, 
#                            alpha, beta, minn,
#                            fdr, thr,
#                            peelcriterion,
#                            seed)
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

cv.tune <- function(traindata, traintime, trainstatus, 
                    alpha, beta, minn, 
                    fdr, thr,
                    peelcriterion,
                    seed) {

  if (!is.null(seed))
    set.seed(seed)

  # Constants
  n <- nrow(traindata)                                   # Number of samples
  p <- ncol(traindata)                                   # Number of covariates
  beta <- max(minn/n, beta)                              # Minimal box support

  # Initialization of vector of pre-selected covariates
  varsel <- 1:p
  names(varsel) <- colnames(traindata)

  # Pre-selection of covariates by univariate SBH
  pval <- numeric(p)
  stat <- numeric(p)
  univarsign <- numeric(p)
  for (j in 1:p) {
    # cat("covariate: ", j, "\r\a")
    prsp.fit <- prsp(traindata=traindata[, j, drop=FALSE], 
                     traintime=traintime, 
                     trainstatus=trainstatus,
                     L=NULL, 
                     varsign=NULL, 
                     alpha=alpha, 
                     beta=beta,
                     peelcriterion=peelcriterion)
    univarsign[j] <- prsp.fit$varsign
    pval[j] <- prsp.fit$pval
    if (peelcriterion == "lr") {
        stat[j] <- prsp.fit$lrt
    } else if (peelcriterion == "hr") {
        stat[j] <- prsp.fit$lhr
    } else {
        stop("Invalid peeling criterion. Exiting...\n")
    }
  }
  if (is.null(thr)) {
    qval <- p.adjust(p=pval, method="BH")
    ord <- order(pval, decreasing=FALSE, na.last=TRUE)
    pval <- pval[ord]
    qval <- qval[ord]
    w <- which(qval <= fdr) 
    if (is.empty(w)) {
        maxsteps <- 0
        varsel <- NULL
        varsign <- NULL
        boxcut <- NULL
    } else {
        varsel <- varsel[ord[w]]
        traindata <- traindata[, varsel, drop=FALSE]
        model.fit <- prsp(traindata=traindata, 
                          traintime=traintime, 
                          trainstatus=trainstatus, 
                          L=NULL, 
                          varsign=univarsign[varsel], 
                          alpha=alpha,
                          beta=beta, 
                          peelcriterion=peelcriterion)
        maxsteps <- model.fit$maxsteps
        varsign <- model.fit$varsign
        boxcut <- model.fit$boxcut
    }
  } else if (is.null(fdr)) {
    ord <- order(stat, decreasing=TRUE, na.last=TRUE)   
    w <- 1:min(p,thr)
    if (thr <= 0) {
        maxsteps <- 0
        varsel <- NULL
        varsign <- NULL
        boxcut <- NULL
    } else {
        varsel <- varsel[ord[w]]
        traindata <- traindata[, varsel, drop=FALSE]
        model.fit <- prsp(traindata=traindata, 
                          traintime=traintime, 
                          trainstatus=trainstatus, 
                          L=NULL, 
                          varsign=univarsign[varsel], 
                          alpha=alpha, 
                          beta=beta, 
                          peelcriterion=peelcriterion)
        maxsteps <- model.fit$maxsteps
        varsign <- model.fit$varsign
        boxcut <- model.fit$boxcut
    }
  }

  # Returning the final results
  return(list("maxsteps"=maxsteps,
              "varsel"=varsel,
              "varsign"=varsign,
              "boxcut"=boxcut))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    prsp (traindata, traintime, trainstatus, L, varsign, alpha, beta, peelcriterion)
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

prsp <- function(traindata, traintime, trainstatus, L, varsign, alpha, beta, peelcriterion) {

  n <- nrow(traindata)                                   # Number of training samples
  p <- ncol(traindata)                                   # Number of covariates
  digits <- getOption("digits")
  ncut <- ceiling(log(beta) / log(1 - (1/n)))            # Maximal number of peeling steps
    
  # Directions of peeling
  if (is.null(varsign)) {
    varsign <- numeric(p)
    for (j in 1:p) {
        # cat("covariate: ", j, "\r\a")
        varsign[j] <- sign(coxph(Surv(traintime, trainstatus) ~ traindata[, j, drop=FALSE], eps=0.01, singular.ok=T, iter.max=1)$coef)
    }  
  }
  names(varsign) <- colnames(traindata)
  
  # Initializations of box boundaries
  initcutpts <- numeric(p)
  for (j in 1:p) {
    if (varsign[j] == 1) {
        initcutpts[j] <- min(traindata[,j])
    } else if (varsign[j] == -1) {
        initcutpts[j] <- max(traindata[,j])
    } else {
        cat("Note: direction of peeling for covariate: ", j, " : ", varsign[j], "\n", sep="")
    }
  }

  # Initializations of returned objects
  vartrace <- numeric(ncut)
  boxcut <- matrix(data=NA, nrow=ncut, ncol=p)
  pvalj <- matrix(NA, ncut, p)
  lhrlj <- matrix(NA, ncut, p)
  lrtlj <- matrix(NA, ncut, p)
  
  # Initializations for the PRSP loop
  boxcutpts <- initcutpts
  boxmass <- 1
  boxes <- matrix(data=FALSE, nrow=n, ncol=p)            # Initial logical matrix of box membership indicator by dimension
  sel <- rep(TRUE, n)                                    # Initial logical vector of in-box samples
  xsel <- traindata                                      # Initial selection of samples from training data
  varpeel <- (apply(traindata, 2, "var") > 10^(-digits)) # Initial selection of covariates for peeling
  continue <- TRUE
  if (!(is.null(L))) {
    switch <- 1
  } else {
    L <- 1
    switch <- 0
  }
  l <- 0

  # PRSP loop
  while ((boxmass >= beta) & (l*switch < L) & (continue)) {
  
    l <- l + 1
    xsign <- t(t(xsel) * varsign)

    # Potential cutpts by dimension
    cutpts.sign <- updatecut(x=xsign, fract=alpha)
    cutpts <- cutpts.sign * varsign

    # Update box membership indicator by dimension
    boxes <- as.matrix(t((t(traindata) * varsign) >= as.vector(cutpts.sign)) & sel)

    vmd <- rep(NA, p)
    for (j in 1:p) {
      boxes1j <- 1 * boxes[,j]
      if ((sum(boxes1j) != length(boxes1j)) && (sum(boxes1j) != 0)) {
        # Rate of increase of LHR (between in and out box)
        if (peelcriterion == "hr") {
          fit <- coxph(Surv(traintime, trainstatus) ~ 1 + boxes1j, singular.ok=TRUE, iter.max=1)
          lhrlj[l,j] <- fit$coef
          pvalj[l,j] <- summary(fit)$coef[,"Pr(>|z|)"]
          if (l == 1) {
            vmd[j] <- (lhrlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lhrlj[l,j] - lhrlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
          # Rate of increase of LRT (between in and out box)
        } else if (peelcriterion == "lr") {
          fit <- survdiff(Surv(traintime, trainstatus) ~ 1 + boxes1j, rho=0)
          lrtlj[l,j] <- fit$chisq
          pvalj[l,j] <- 1 - pchisq(q=fit$chisq, df=1)
          if (l == 1) {
            vmd[j] <- (lrtlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lrtlj[l,j] - lrtlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        } else {
          stop("Invalid peeling criterion. Exiting...\n")
        }
      } else {
        varpeel[j] <- FALSE
      }
    }

    # If the previous attempted peeling succeeded
    if (sum(varpeel) > 0 && (!is.empty(vmd[(!is.nan(vmd)) & (!is.infinite(vmd)) & (!is.na(vmd))]))) {
      
      # Maximizing the rate of increase of LHR or LRT (peeling criterion).
      # Only one variable (the first one in rank) is selected in case of ties
      varj <- which(vmd == max(vmd[(!is.nan(vmd)) & (!is.infinite(vmd))], na.rm=TRUE))[1]

      # Updating
      sel <- boxes[, varj]
      boxmass <- mean(1 * sel)
      xsel <- traindata[sel, ,drop=FALSE]
      varpeel <- (apply(xsel, 2, "var") > 10^(-digits))
      boxcutpts[varj] <- cutpts[varj]
      
      # Saving trained box quantities of interest for the current peeling step
      boxcut[l, ] <- boxcutpts
      vartrace[l] <- varj
      
    # Else exiting the loop and decrementing the peeling step number since the last attempted failed in that case
    } else {
      continue <- FALSE
      l <- l - 1
    }
    
  }
  pval <- pvalj[l,,drop=TRUE]
  lrt <- lrtlj[l,,drop=TRUE]
  lhr <- lhrlj[l,,drop=TRUE]

  if (l == 0) {
    # Taking the first step box covering all the data
    boxcut <- as.matrix(initcutpts)
    vartrace <- 0
  } else if (l >= 1) {
    # Prepending the first step box covering all the data
    boxcut <- rbind(initcutpts, boxcut[1:l, , drop=FALSE])
    vartrace <- c(0, vartrace[1:l])
  }
  rownames(boxcut) <- paste("step", 0:l, sep="")
  colnames(boxcut) <- colnames(traindata)
  names(vartrace) <- paste("step", 0:l, sep="")
  l <- l + 1  

  return(list("maxsteps"=l, "boxcut"=boxcut, "vartrace"=vartrace, "varsign"=varsign, "lrt"=lrt, "lhr"=lhr, "pval"=pval))
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    cv.folds (n, K, seed=NULL)
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

cv.folds <- function (n, K, seed=NULL) {

  if (!is.null(seed))
    set.seed(seed)

  n <- round(rep(n, length.out = 1))
  if (!isTRUE(n > 0))
    stop("'n' must be positive")
  K <- round(rep(K, length.out = 1))
  if (!isTRUE((K >= 1) && K <= n))
    stop(paste("'K' outside allowable range {1,...,", n, "} \n", sep=""))
  if (K == 1) {
    perm <- seq_len(n)
  } else if (K == n) {
    perm <- seq_len(n)
  } else {
    perm <- sample(n)
  }
  permkey <- pmatch(x=1:n, table=perm)
  w <- rep(seq_len(K), length.out=n)
  ord <- numeric(0)
  for (k in 1:K) {
    ord <- c(ord, perm[(w == k)])
  }
  foldkey <- pmatch(x=1:n, table=ord)
  folds <- list(n=n, K=K, perm=perm, permkey=permkey, which=w, foldkey=foldkey, seed=seed)

  return(folds)
}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    endpoints (ind, timemat, probmat, timeval, probval)
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

endpoints <- function(ind, timemat, probmat, timeval, probval) {

  N <- nrow(timemat) #N <- nrow(probmat)
  L <- N-length(ind)
  if (!(is.empty(ind))) {
    timemat <- timemat[-ind, , drop=FALSE]
    probmat <- probmat[-ind, , drop=FALSE]
  }
  min.prob.bar <- apply(probmat, 1, min, na.rm=TRUE)
  max.time.bar <- apply(timemat, 1, max, na.rm=TRUE)
  if (is.null(probval) && is.null(timeval)) {
    prob.bar <- rep(NA, L)
    time.bar <- rep(NA, L)
  } else if (!is.null(probval)) {
    prob.bar <- rep(probval, L)
    ind.probmat <- (probmat <= probval)
    ind.probmat[is.na(ind.probmat)] <- FALSE
    time.bar <- numeric(L)
    for (l in 1:L) {
      if (probval >= min.prob.bar[l]) {
        time.bar[l] <- min(timemat[l,which(ind.probmat[l,,drop=TRUE])])
      } else {
        time.bar[l] <- NA
      }
    }
  } else if (!is.null(timeval)) {
    time.bar <- rep(timeval, L)
    ind.timemat <- (timemat >= timeval)
    ind.timemat[is.na(ind.timemat)] <- FALSE
    prob.bar <- numeric(L)
    for (l in 1:L) {
      if (timeval <= max.time.bar[l]) {
        prob.bar[l] <- max(probmat[l,which(ind.timemat[l,,drop=TRUE])])
      } else {
        prob.bar[l] <- NA
      }
    }
  }

  return(list("time.bar"=time.bar,
              "prob.bar"=prob.bar,
              "max.time.bar"=max.time.bar,
              "min.prob.bar"=min.prob.bar))

}
##########################################################################################################################################




##########################################################################################################################################
#################
# Usage         :
################
#                    updatecut (x, fract)
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

updatecut <- function(x, fract) {

  p <- dim(x)[2]
  cutpts <- apply(x, 2, "quantile", type=7, probs=fract)

  for (j in 1:p) {
    xunique <- sort(unique(x[, j]))
    if (length(xunique) == 1) {
      cutpts[j] <- min(xunique)
    } else {
      if (isTRUE(all.equal(as.single(cutpts[j]), as.single(min(xunique))))) {
        cutpts[j] <- xunique[2]
      }
    }
  }

  return(cutpts)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   lapply.array (X, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA, MARGIN=1:2, FUN, ...)
#
################
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

lapply.array <- function (X, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA, MARGIN=1:2, FUN, ...) {
  x <- list2array(list=X, rowtrunc=rowtrunc, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   lapply.mat (X, coltrunc=NULL, sub=NULL, fill=NA, MARGIN=2, FUN, ...)
#
################
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

lapply.mat <- function (X, coltrunc=NULL, sub=NULL, fill=NA, MARGIN=2, FUN, ...) {
  x <- list2mat(list=X, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   list2array (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA)
#
################
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

list2array <- function (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    min.row <- min(sapply(my.list, nrow))
    max.row <- max(sapply(my.list, nrow))
    min.col <- min(sapply(my.list, ncol))
    max.col <- max(sapply(my.list, ncol))
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[,1:min.col,drop=FALSE]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
      } else {
        if (coltrunc <= min.col) {
            adjusted.list <- lapply(my.list, function(x) {x[,1:coltrunc,drop=FALSE]})
        } else if ((coltrunc > min.col) & (coltrunc <= max.col)) {
            adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max(0,coltrunc - ncol(x))))})
        } else {
            adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=coltrunc - ncol(x)))})
        }
      } 
    } else {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
    }
    if (!is.null(rowtrunc)) {
      if (rowtrunc == "min") {
        adjusted.list <- lapply(adjusted.list, function(x) {x[1:min.row,,drop=FALSE]})
      } else if (rowtrunc == "max") {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
      } else {
        if (rowtrunc <= min.row) {
            adjusted.list <- lapply(adjusted.list, function(x) {x[1:rowtrunc,,drop=FALSE]})
        } else if ((rowtrunc > min.row) & (rowtrunc <= max.row)) {
            adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max(0,rowtrunc - nrow(x)), ncol=ncol(x)))})
        } else {
            adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=rowtrunc - nrow(x), ncol=ncol(x)))})
        }
      } 
    } else {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
    }
    my.array <- array(data=fill, dim=c(nrow(adjusted.list[[1]]), ncol(adjusted.list[[1]]), length(adjusted.list)))
    for(i in 1:length(adjusted.list)) {
      my.array[,,i] <- adjusted.list[[i]]
    }
  } else {
    my.array <- array(data=fill, dim=c(0,0,0))
  }
  return(my.array)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   list2mat (list, coltrunc=NULL, sub=NULL, fill=NA)
#
################
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

list2mat <- function (list, coltrunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    min.col <- min(sapply(my.list, length))
    max.col <- max(sapply(my.list, length))
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[1:min.col]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, max.col - length(x)))})
      } else {
        if (coltrunc <= min.col) {
            adjusted.list <- lapply(my.list, function(x) {x[1:coltrunc]})
        } else if ((coltrunc > min.col) & (coltrunc <= max.col)) {
            adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, max(0,coltrunc - length(x))))})
        } else {
            adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, coltrunc - length(x)))})
        }
      } 
    } else {
        adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, max.col - length(x)))})
    }
    my.mat <- do.call(rbind, adjusted.list)
  } else {
    my.mat <- matrix(data=fill, nrow=0, ncol=0)
  }
  return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   cbindlist (list, trunc)
#
################
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

cbindlist <- function(list, trunc) {
  if (!is.empty(list)) {
    max.row <- max(sapply(list, nrow))
    adjusted.list <- lapply(list, function(x) {rbind(x, matrix(data=NA, nrow=max.row - nrow(x), ncol=ncol(x)))})
    my.mat <- adjusted.list[[1]]
    lcl <- length(adjusted.list)
    if (lcl > 1) {
      for(i in 2:lcl){
        my.mat <- cbind(my.mat, adjusted.list[[i]])
      }
    }
    if (missing(trunc)) {
      trunc <- max.row
    }
    my.mat <- my.mat[1:trunc,,drop=FALSE]
  } else {
    my.mat <- matrix(data=NA, nrow=0, ncol=0)
  }
  return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   is.empty(x)
################
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

is.empty <- function(x) {
  if (is.vector(x)) {
    if((length(x) == 0) || (x == "")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (is.matrix(x) || is.data.frame(x)) {
    return( ((nrow(x) == 0) || (ncol(x) == 0)) )
  } else {
    return( ((length(x) == 0) || (x == "")) )
  }
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   myround (x, digits = 0)
################
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

myround <- function (x, digits = 0) {
    upround <- function (x, digits = 0) {
        return(ceiling(x*10^(digits))/10^digits)
    }
    dnround <- function (x, digits = 0) {
        return(floor(x*10^digits)/10^digits)
    }
    i <- (x - trunc(x) >= 0.5)
    x[i] <- upround(x[i], digits=digits)
    x[!i] <- dnround(x[!i], digits=digits)
    return(x)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :   .onAttach (libname, pkgname)
################
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

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("Type PRIMsrc.news() to see the log file of new features, changes, and bug fixes\n")
}
##########################################################################################################################################
