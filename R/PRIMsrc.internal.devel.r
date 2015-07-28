##########################################################################################################################################
# PRIMsrc
##########################################################################################################################################

##########################################################################################################################################
# SURVIVAL INTERNAL SUBROUTINES
# (never to be called by end-user)
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
