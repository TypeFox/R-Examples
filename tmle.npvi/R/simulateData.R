simulateData <- function(B, W, X, Xq, g, mu, muAux, sigma2, theta=NULL, Y=list(value=NA, index=NA), weightsW=rep(1, length(W)), family=c("parsimonious", "gaussian"), verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'B':
  B <- Arguments$getNumeric(B);
  
  ## Argument 'W':
  W <- Arguments$getNumerics(W);

  ## Argument 'X':
  X <- Arguments$getNumerics(X);

  ## Argument 'Xq':
  Xq.value <- Arguments$getNumerics(Xq$value);
  Xq.index <- Arguments$getIntegers(Xq$index);
  nMax <- nrow(Xq)
  
  ## Argument 'g':
  mode <- mode(g);
  if (mode != "function") {
    throw("Argument 'g' should be of mode 'function', not '", mode);
  }

  ## Argument 'mu':
  mode <- mode(mu);
  if (mode != "function") {
    throw("Argument 'mu' should be of mode 'function', not '", mode);
  }

  ## Argument 'muAux':
  mode <- mode(muAux);
  if (mode != "function") {
    throw("Argument 'muAux' should be of mode 'function', not '", mode);
  }

  
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);
  if (is.na(sigma2)) {
    throw("Argument 'sigma2' should be provided...")
  }

  ## Argument 'theta'
  if (!is.null(theta)) {
    mode <- mode(theta)
    if (mode != "function") {
      throw("Argument 'theta' should be 'NULL' or of mode 'function', not '",  mode)
    }
  }
  
  ## Argument 'family':
  family <- match.arg(family);

  ## Argument 'Y'
  Y.value <- Arguments$getNumerics(Y$value);
  Y.index <- Arguments$getNumerics(Y$index);
  if (!is.null(theta)) {
    if (any(is.na(Y.value)) & family=="gaussian") {
      throw("Argument 'Y$value' of mode 'numerics' should be provided when 'family' is 'gaussian'.")
    }
  }
    
  ## Argument 'weightsW':
  weightsW <- Arguments$getNumerics(weightsW);
  nr <- length(W)
  if (length(weightsW) != nr) {
    throw("Length of 'weightsW' must match length of W!")
  }

  ## Argument 'verbose':
  verbose <- Arguments$getLogical(verbose);

  ## the probability P(X=0) directly computed on the whole dataset
  meanGW <- mean(X==0);
  whichXisZero <- which(X==0);
  whichXisNotZero <- which(X!=0);
  obsX <- X[whichXisNotZero];
  rm(X);

  WB <- sample(W, B, replace=TRUE, prob=weightsW)
  XB <- rep(NA, B)
  YB <- rep(NA, B)
  muWB <- mu(WB)
  muAuxWB <- muAux(WB)
  gWB <- g(WB)
  ##
  U <- (runif(B) >= gWB)
  if (family=="gaussian") {
    XB[!U] <- 0
  } else if (family=="parsimonious") {
    ## old:
    ## XB[!U] <- whichXisZero[1] ## first index of row with X equal to 0
    ## new:
    XB[!U] <- Xq.index[Xq.value==0]
  }
  ##
  gW <- gWB[U]
  condMeanX <- muAuxWB[U]

  ##
  parameters <- list(meanGW=meanGW, muWB=muWB, gWB=gWB, U=U)
  if (family=="gaussian") {
    ## Note: here, 'tabulate' is necessarily FALSE
    ## hence 'W' are actual observations and not indices
    if (FALSE) {
      sigma2Bis <- sigma2 - mean(muWB^2/(1-gWB))
    } else {
      sigma2Bis <- sigma2 - mean(muAuxWB^2*(1-gWB))
    }
    if (sigma2Bis <=0) {
      cat("sigma2:\n")
      print(sigma2)
      cat("mean(muAuxWB^2*(1-gWB)):\n")
      print(mean(muAuxWB^2*(1-gWB)))
      throw("Problem here!... Parameter 'sigma2Bis' should be positive...")
    }
    condVarX <- sigma2Bis/(1-gW)
    XB[U] <- rnorm(sum(U), mean=condMeanX, sd=sqrt(condVarX))
    if (!is.null(theta)) {
      YB <- rnorm(B, mean=theta(cbind(X=XB, W=WB)), sd=sd(Y.value))
    }
  } else if (family=="parsimonious") {
    indices <- simulateParsimoniouslyXgivenW(WB[U], min(obsX), max(obsX),
                                             Xq, condMeanX, sigma2, parameters)
    ## ## CAUTION
    Xq.index <- Xq.index[Xq.value!=0]
    XB[U] <- Xq.index[indices]
    if (!is.null(theta)) {
      T <- theta(cbind(X=XB, W=WB))
      indices <- simulateParsimoniouslyYgivenXW(T, Y.value)
      YB <- Y.index[indices]
    }
  }
  obsB <- cbind(W=WB, X=XB, Y=YB)
  rm(WB, XB, YB, U)
  obsB
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

