######################################################################
# Compute log likelihood ratio for a univariate or multivariate
# categorical distribution
#
# Params:
#  outcomes - a data frame with all possible configuration for the (c-1)
#             variables not being the reference category.
#  mu  - expectation under which LLR under pi is computed
#  mu0 - null model. A vector of length (k-1)
#  mu1 - alternative model. A vector of length (k-1)
######################################################################

LLR.fun <- function(outcomes, mu, mu0, mu1, dfun, ...) {  
  #Compute likelihood ratios. Both univariate and the multivariate
  #values are computed
  llr.res <-   t(apply(outcomes,1, function(y) {
    llr <- dfun(y, mu=mu1, log=TRUE,...) - dfun(y, mu=mu0, log=TRUE, ...)
    p <- dfun(y, mu=mu, ...)
   return(c(llr=llr,p=p)) 
  }))
  res <- cbind(outcomes,llr.res)
  colnames(res) <- c(paste("y",1:ncol(outcomes),sep=""),"llr","p")
  return(res)
}

######################################################################
# Function to compute all possible outcomes for the categorical time
# series. This is needed for the LLR computations
#
# Parameters:
# km1 - Dimension of the problem (k-1)
# n   - number of items arranged (i.e. number of experiments). Integer
#
# Returns:
#  matrix of size (number of configs) \times km1
#  containing all possible states
######################################################################

outcomeFunStandard <- function(k,n) {
  #Compute all possible likelihood ratios and their probability under mu
  #Note: Currently all states are investigated. This might be way too
  #much work as defacto many states have an occurence prob near 0!!
  args <- list() ; for (j in seq_len(k)) args[[j]] <- 0:n
  outcomes <- as.matrix(do.call("expand.grid", args))
  
  #Take only valid outcomes (might reduce drastically the number of cells)
  if (!is.null(n)) {
    outcomes <- outcomes[apply(outcomes,1,sum) <= n,,drop=FALSE]
  }
  return(outcomes)
}

######################################################################
# Compute run length for CUSUM based on Markov representation of the
# Likelihood ratio based CUSUM
#
# Parameters:
#  mu -  (k-1 \times T) matrix with true proportions, i.e. equal to mu0 or mu1 if one wants to compute e.g. ARL_0 or ARL_1
#  mu0 - (k-1 \times T) matrix with in-control proportions
#  mu1 - (k-1 \times T) matrix with out-of-control proportion 
#  n   - vector of length T containing the total number of experiments for each time point
#  h- The threshold h which is used for the CUSUM
#  g   - The number of levels to cut the state space into, i.e. M on foil 12
######################################################################

LRCUSUM.runlength <- function(mu,mu0,mu1,h,dfun, n, g=5,outcomeFun=NULL,...) {
  #Semantic checks
  if ( ((ncol(mu) != ncol(mu0)) | (ncol(mu0) != ncol(mu1))) |
      ((nrow(mu) != nrow(mu0)) | (nrow(mu0) != nrow(mu1)))) {
    stop("Error: dimensions of mu, mu0 and mu1 have to match")
  }
  if (missing(h)) {
    stop("No threshold specified!")
  }
  #If no specific way for computing the outcomes is given
  #use the standard way.
  if (is.null(outcomeFun)) {
    outcomeFun <- outcomeFunStandard
  }
  
  #Discretize number of possible states of the CUSUM
  S <- c(-Inf,seq(0,h,length=g))
  names <- c(levels(cut(1,S,right=TRUE)),">=h")
  #Time variable
  t <- 1:ncol(mu)
  #Dimension of the problem (k-1)
  km1 <- nrow(mu)
  
  #Create transition matrix for CUSUM control chart
  P <- array(0, dim=c(length(t),g+1,g+1),dimnames=list(t,names,names))
  #Once in the absorbing state stay there!
  P[,g+1,g+1] <- 1
  
  #Loop over all P[t,,] and compute probabilities
  for (i in seq_len(length(t))) {
    cat("Looking at t=",i," out of ",length(t),"\n")

    #Determine all possible outcomes
    outcomes <- outcomeFun(km1,n[i])

    #Compute all possible likelihood ratios and their probability under mu
    llr <- LLR.fun(outcomes,mu=mu[,i],mu0=mu0[,i],mu1=mu1[,i],dfun=dfun,size=n[i],...)

    #Exact CDF of the LLR for this time
    F <- stepfun(sort(llr[,"llr"]),c(0,cumsum(llr[order(llr[,"llr"]),"p"])))

    #Compute probability going from c <= S_{t-1} < d to a <= S_{t} < b
    for (j in 1:g) { #from index
      for (k in 1:g) { #to index
        a <- S[k] ; b <- S[k+1] ; c <- S[j] ; d <- S[j+1] ; m <- (c+d)/2

        #From zero to new state
        if (j == 1) {
          P[i,j,k] <- F(b) - F(a)
        } else { 
            #Rieman integral assuming as in Brook & Evans (1972) that S at midpoint
            #P[i,j,k] <- F(b-m) - F(a-m)
            #Slightly better approximation by Hawkins (1992), which uses Simpson's rule
            P[i,j,k] <- (F(b-c) + 4*F(b-m) + F(b-d) - F(a-c) - 4*F(a-m) - F(a-d))/6
        }
      }
 }
    #Whatever is missing goes to >h category (take care of rounding errors)
    P[i,-(g+1),(g+1)] <- pmax(0,1-apply(P[i,-(g+1),-(g+1)],1,sum))
  }

  #Use matrix to compute RL distribution
  Ppower <- P[1,,]
  alarmUntilTime <- numeric(ncol(mu0))
  alarmUntilTime[1] <- Ppower[1,ncol(P)]
  for (time in t[-1]) { #from 2 to length of t
    Ppower <- Ppower %*% P[time,,]
    alarmUntilTime[time] <- Ppower[1,ncol(P)]
  }
  pRL <- c(alarmUntilTime[1],diff(alarmUntilTime))

  mom <- NA
  #If the Markov chain is homogenous then compute ARL by inverting
  if (length(t) == 1) {
    R <- P[,1:g,1:g]
    I <- diag(rep(1,g))
    mom <- solve(I-R) %*% matrix(1,g,1) #-- no sparse computing
    #Alternative using sparse 
    #mom <- solve(Matrix(id-r)) %*% matrix(1,mw,1)
  }
       
  return(list(P=P,pmf=pRL,cdf=alarmUntilTime,arl=mom[1]))
}
