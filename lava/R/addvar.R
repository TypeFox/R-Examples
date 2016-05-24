##' Generic method for adding variables to model object
##'
##' @title Add variable to (model) object
##' @param x Model object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
##' @aliases addvar<-
##' @export
`addvar` <-
function(x,...) UseMethod("addvar")

##' @export
`addvar<-` <-
function(x,...,value) UseMethod("addvar<-")


##' @export
`addvar<-.lvm` <- function(x,...,value) {
  if (inherits(value,"formula")) {
    regression(x,...) <- value
    return(x)
##    return(addvar(x,all.vars(value),...))
  }
  addvar(x, var=value, ...)
}

##' @export
`addvar.lvm` <-
function(x, var, silent=lava.options()$silent,reindex=TRUE,...) {
  new <- setdiff(var,vars(x))
  k <- length(new)
  Debug(new)
  if (k>0) {
    if (lava.options()$sparse) {
      requireNamespace("Matrix",quietly=TRUE)
      newNA <- newM <- Matrix::Matrix(0,k,k)
      newNAc <- newNA; diag(newNAc) <- NA
      newcov <- Matrix::Diagonal(k)
    } else {
      newM <- matrix(0,k,k)
      newcov <- diag(k)
    }
    newNA <- matrix(NA,k,k)
    colnames(newM) <- rownames(newM) <-
      colnames(newcov) <- rownames(newcov) <-
        colnames(newNA) <- rownames(newNA) <- new
    newmean <- as.list(rep(NA,k))
    ##  for (i in new) {
    N <- nrow(x$cov)
    if (is.null(N)) {
      N <- 0
      x$M <- newM
      x$cov <- newcov; x$covfix <- x$fix <- x$par <- x$covpar <- newNA
      x$mean <- newmean
    } else {
      if (lava.options()$sparse) {
        x$M <- Matrix::bdiag(x$M, newM) ## Add regression labels.R
        x$cov <- Matrix::bdiag(x$cov, newcov) ## Add covariance
        x$par <- Matrix::bdiag(x$par, newNA) ## Add regression labels
        x$covpar <- Matrix::bdiag(x$covpar, newNA) ## Add covariance labels
        x$fix <- Matrix::bdiag(x$fix, newNA)
        x$covfix <- Matrix::bdiag(x$covfix,  newNA)

      } else {
        x$M <- blockdiag(x$M, newM, pad=0) ## Add regression labels
        x$cov <- blockdiag(x$cov, newcov, pad=0) ## Add covariance
        x$par <- blockdiag(x$par, newNA, pad=NA) ## Add regression labels
        x$covpar <- blockdiag(x$covpar, newNA, pad=NA) ## Add covariance labels
        x$fix <- blockdiag(x$fix, newNA, pad=NA) ##
        x$covfix <- blockdiag(x$covfix,  newNA, pad=NA) ##
      }
      x$mean <- c(x$mean, newmean)
    }
    names(x$mean)[N+seq_len(k)] <-
      colnames(x$M)[N+seq_len(k)] <- rownames(x$M)[N+seq_len(k)] <-
        colnames(x$covfix)[N+seq_len(k)] <- rownames(x$covfix)[N+seq_len(k)] <-
          colnames(x$fix)[N+seq_len(k)] <- rownames(x$fix)[N+seq_len(k)] <-
            colnames(x$covpar)[N+seq_len(k)] <- rownames(x$covpar)[N+seq_len(k)] <-
              colnames(x$par)[N+seq_len(k)] <- rownames(x$par)[N+seq_len(k)] <-
                colnames(x$cov)[N+seq_len(k)] <- rownames(x$cov)[N+seq_len(k)] <- new

    ## x$cov[N+1,N+1] <- 1
    ## names(x$mean)[N+1] <-
    ##   colnames(x$M)[N+1] <- rownames(x$M)[N+1] <-
    ##     colnames(x$covfix)[N+1] <- rownames(x$covfix)[N+1] <-
    ##       colnames(x$fix)[N+1] <- rownames(x$fix)[N+1] <-
    ##         colnames(x$covpar)[N+1] <- rownames(x$covpar)[N+1] <-
    ##           colnames(x$par)[N+1] <- rownames(x$par)[N+1] <-
    ##             colnames(x$cov)[N+1] <- rownames(x$cov)[N+1] <- i
    ## myexpr <- paste("c(",i,"=expression(",i,"))", sep="\"")
    ## labels(x) <- (eval(parse(text=myexpr)))
    ## if (!silent)
    ##   message("\tAdded '", i, "' to model.\n", sep="")
    if (!silent) {
      if (k==1)
        message("\tAdded '", new, "' to model.\n", sep="")
      else
        message("\tAdded ",paste(paste("'",new,"'",sep=""),collapse=",")," to model.\n", sep="")
    }
    exogenous(x) <- c(new,exogenous(x))
  }
  if (reindex)
    index(x) <- reindex(x)
  return(x)
}
