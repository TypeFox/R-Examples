##'Composite Likelihood for probit latent variable models
##'
##'Estimate parameters in a probit latent variable model via a composite
##'likelihood decomposition.
##'
##'
##'@param x \code{lvm}-object
##'@param data data.frame
##'@param k Size of composite groups
##'@param type Determines number of groups. With \code{type="nearest"} (default)
##'only neighboring items will be grouped, e.g. for \code{k=2}
##'(y1,y2),(y2,y3),... With \code{type="all"} all combinations of size \code{k}
##'are included
##'@param pairlist A list of indices specifying the composite groups. Optional
##'argument which overrides \code{k} and \code{type} but gives complete
##'flexibility in the specification of the composite likelihood
##'@param silent Turn output messsages on/off
##'@param \dots Additional arguments parsed on to lower-level functions
##'@return An object of class \code{clprobit} inheriting methods from \code{lvm}
##'@author Klaus K. Holst
##'@seealso \link[lava]{estimate}
##'@keywords models regression
##' @export
clprobit <- function(x,data,k=2,type=c("nearest","all"),pairlist,silent=TRUE,
                     ...) {
  y <- endogenous(x)
  binsurv <- rep(FALSE,length(y))
  for (i in 1:length(y)) {
    z <- data[,y[i]]
    binsurv[i] <- is.Surv(z) | (is.factor(z) && length(levels(z))==2)
  }
  
  binsurv <- unique(c(y[binsurv],binary(x)))
  ##  binsurvpos <- which(colnames(data)%in%binsurv)
  if (!missing(pairlist)) {
    binsurvpos <- which(colnames(data)%in%endogenous(x))
  } else {
    binsurvpos <- which(colnames(data)%in%binsurv)
  }
  
  if (missing(pairlist)) {
    if (length(binsurv)<(k+1)) stop("No need for composite likelihood analysis.")

    if (type[1]=="all") {
      mypar <- combn(length(binsurv),k) ## all pairs (or multiplets), k=2: k*(k-1)/2
    } else {
      mypar <- sapply(0:(length(binsurv)-k), function(x) x+1:k)
    }
  } else {
    mypar <- pairlist
  }  
  
  if (is.matrix(mypar)) {
    mypar0 <- mypar; mypar <- c()
    for (i in seq(ncol(mypar0)))
      mypar <- c(mypar, list(mypar0[,i]))
  }
  
  nblocks <- length(mypar)
  mydata0 <- data[c(),,drop=FALSE]  
  mydata <-  as.data.frame(matrix(NA, nblocks*nrow(data), ncol=ncol(data)))
  names(mydata) <- names(mydata0)
  for (i in 1:ncol(mydata)) {
    if (is.factor(data[,i])) {
      mydata[,i] <- factor(mydata[,i],levels=levels(mydata0[,i]))
    }
    if (is.Surv(data[,i])) {
      S <- data[,i]
      for (j in 2:nblocks) S <- rbind(S,data[,i])
      S[,1] <- NA
      mydata[,i] <- S
    }
  }
  
  for (ii in 1:nblocks) {    
    data0 <- data;
    for (i in binsurvpos[-mypar[[ii]]]) {
      if (is.Surv(data[,i])) {
        S <- data0[,i]; S[,1] <- NA
        data0[,i] <- S
      } else {
        data0[,i] <- NA
        if (is.factor(data[,i])) data0[,i] <- factor(data0[,i],levels=levels(data[,i]))
      }
    }
    mydata[(1:nrow(data))+(ii-1)*nrow(data),] <- data0
##    mydata <- rbind(mydata,data0)
  }

  suppressWarnings(e0 <- estimate(x,data=mydata,missing=TRUE,silent=silent,
              ...))

  S <- score(e0,indiv=TRUE)
  nd <- nrow(data)
  block1 <- which((1:nd)%in%(rownames(S)))
  blocks <- sapply(1:nblocks, function(x) 1:length(block1)+length(block1)*(x-1))
  Siid <- matrix(0,nrow=length(block1),ncol=ncol(S))
  for (j in 1:ncol(blocks)) {
    Siid <- Siid+S[blocks[,j],]
  }
  iI <- vcov(e0); J <- t(Siid)%*%(Siid)
  e0$iidscore <- Siid
  e0$blocks <- blocks
  e0$vcov <- iI%*%J%*%iI ## thetahat-theta0 :=(asymp) I^-1*S => var(thetahat) = iI*var(S)*iI 
  cc <- e0$coef; cc[,2] <- sqrt(diag(e0$vcov))
  cc[,3] <- cc[,1]/cc[,2]; cc[,4] <- 2*(1-pnorm(abs(cc[,3])))
  e0$coef <- cc
  class(e0) <- c("clprobit",class(e0))
  return(e0)
}
score.clprobit <- function(x,indiv=FALSE,...) {
  if (!indiv)
    return(colSums(x$iidscore))
  x$iidscore
}

