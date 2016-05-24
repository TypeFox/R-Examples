#' @title Quick and dirty stochastic generation of seasonal streamflow replicates for a single site.
#' @description Generates seasonal time series using either the kNN Bootstrap (non-parametric) or a numerically-fitted PARMA(1,1) (parametric) model. For the parametric model, the function automatically transforms the seasonal sub-series to normal and deseasonalizes prior to model fitting.
#' @param Q             time series object with seasonal resolution (e.g., frequency = 2, 3, 4, 6 or 12 for monthly data). 
#' @param reps          integer. The number of replicates to be generated.  The default is 100.
#' @param years         integer. The length of each replicate in years. The default is equal to the number of complete years given in Q.
#' @param k             integer. The k parameter of the kNN Bootstrap (i.e., number of nearest neighbors). If left blank k = n ^ 0.5., where n is the number of years in the input data.
#' @param d             integer. The d parameter of the kNN Bootstrap (i.e., number of previous time periods to inform the model). If left blank d = 1.
#' @param adjust        logical. If TRUE (the default) the final output time series X will be coerced for 0 <= X <= 1.2*max(Q). Applies only if the PARMA method is used.
#' @param parameters    logical. If TRUE the output will be given as a list including the replicate samples and relevant model parameters (k and d for kNNboot and phi, theta and standard deviation of residuals for PARMA). The default is FALSE.
#' @param method        character string giving the method used to generate the data. Defaults to "kNNboot" - the k Nearest Neighbour Bootstrap. See references for detail on the two methods available.
#' @return Returns a multi time series object containing synthetic streamflow replicates.
#' @examples
#' Q <- resX$Q_Mm3
#' replicates <- dirtyreps(Q, reps = 3)
#' mean(replicates); mean(Q)
#' sd(replicates); sd(Q)
#' plot(replicates)
#' @references kNN Bootstrap method: Lall, U. and Sharma, A., 1996. A nearest neighbor bootstrap for resampling hydrologic time series. Water Resources Research, 32(3), pp.679-693.
#' @references PARMA method: Salas, J.D. and Fernandez, B., 1993. Models for data generation in hydrology: univariate techniques. In Stochastic Hydrology and its Use in Water Resources Systems Simulation and Optimization (pp. 47-73). Springer Netherlands.
#' @import stats
#' @importFrom utils head tail
#' @export
dirtyreps <- function(Q, reps, years, k, d, adjust, parameters, method = "kNNboot"){
  
  #if (missing(reps)) stop("reps missing! Enter the required number of replicates. E.g., reps = 100")
  #if (missing(years)) stop("years missing! Enter number of years required for each replicate")

  
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object")
  
  frq <- frequency(Q)
  
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }
  
  len <- length(Q)
  cyc <- cycle(Q)
  
  if (missing(reps)) reps = 100
  if (missing(years)) years = len/frq
  
  if (method == "kNNboot"){
    # CREATE FUNCTION TO GET EUCLIDIAN DISTANCE
    getEucDist <- function(D_t, D_i) {
      EucDist <- sqrt(sum((D_t - D_i)^2))
      return(EucDist)
    }
    # SET UP PARAMETERS
    if (missing(k)) k <- round((len/frq) ^ 0.5)
    if (missing(d)) d <- 1
    # DEFINE PROBABILITY MASS FUNCTION FOR THE KERNAL
    kernal.weights <- (1 / seq(1:k)) / (sum(1 / seq(1:k)))
    # SET UP FOR SIMULATION AND INITIALIZE WITH DATA FROM Q
    nsamples <- years * frq * reps
    X <- vector("numeric", nsamples)
    Xcyc <- rep(1:frq, reps * years)
    X[1:d] <- Q[1:d]
    # SIMULATE THE FULL REPLICATE SEQUENCE
    for (t in (d + 1):length(X)){
      D_i <- X[(t - d):(t - 1)]
      period <- Xcyc[t]
      seqD_t <- Xcyc[(t - d):(t - 1)]
      D_t <- matrix(Q[which(cyc %in% seqD_t)], ncol = d, byrow = TRUE)
      distances <- apply(D_t, 1, getEucDist, D_i = D_i)
      q <- Q[cyc==Xcyc[t]]
      x_pot <- q[match(1:k, rank(distances, ties.method = "random"))]
      X[t] <- sample(x_pot, 1, prob = kernal.weights)
    }
    # SPLIT AND CONVERT TO TIME SERIES
    X <- ts(matrix(X, ncol = reps), frequency = frq)
    params <- c(k, d)
    names(params) <- c("k", "d")
  }
  if (method == "PARMA"){
    # CREATE FUNCTION TO ESTIMATE PARAMETER OF THE LOG TRANSFORM
    getlogtransparam <- function(x){
      a <- (min(x) * max(x) - median(x) ^ 2) / (min(x) + max(x) - 2 * median(x))
      if(a >= min(x)) a <- min(x) - 1
      return(a)
    }
    # CREATE FUNCTION TO LOG-TRANSFORM SERIES
    getQtrans <- function(x){
      a <- (min(x) * max(x) - median(x) ^ 2) / (min(x) + max(x) - 2 * median(x))
      if(a >= min(x)) a <- min(x) - 1
      xT <- log(x - a)
      return(xT)
    }
    # CREATE FUNCTIONS FOR PARMA FITTING
    get_params <- function(x){
      model <- lm(x[,1] ~ 0 + x[,2] + x[,4])
      ph <- model[[1]][[1]] 
      th <- model[[1]][[2]]
      return(c(ph,th))
    }
    get_resids <- function(x){
      model <- lm(x[,1] ~ 0 + x[,2] + x[,4])
      return(model[[2]])
    }
    # TRANSFORM THE SEASONAL SUBSERIES INDIVIDUALLY
    Qtr <- ts(as.vector(t(matrix(unlist(tapply(Q, cyc, getQtrans)),ncol = frq))),
              start = start(Q), frequency = frq)
    ## REMOVE MEAN AND DIVIDE BY STANDARD DEVIATION TO GET WHITE NOISE
    Qtr_ds <- ((Qtr - rep(tapply(Qtr, cyc, mean),len / frq)) / 
                 rep(tapply(Qtr, cyc, sd), len / frq))
    ## INITIALISE RESIDUALS FOR LEAST SQUARES PARMA FITTING
    M <- list()
    Y <- Qtr_ds[cyc == 1][2:(len / frq)]
    X <- Qtr_ds[cyc == frq][1:((len / frq) - 1)]
    ar1_model <- lm(Y ~ 0 + X)
    E <- as.numeric(ar1_model[[2]])
    M[[1]] <- cbind(Y, X, E, vector("numeric", length(X)))
    for(v in 2:frq){
      Y <- Qtr_ds[cyc == v]
      X <- Qtr_ds[cyc == (v - 1)]
      ar1_model <- lm(Y ~ 0 + X)
      E <- as.numeric(ar1_model[[2]])
      M[[v]] <- cbind(Y, X, E, vector("numeric", length(X)))
    }
    paramsx <- matrix(0, frq, 2)
    #ITERATE
    n = 0
    repeat{
      n = n + 1
      for(i in 1:frq){
        if(i==1){
          M[[i]][,4] <- head(M[[frq]][,3], -1)
        }else if(i==2){
          M[[i]][,4] <- c(0 , M[[i-1]][,3])  
        }else{
          M[[i]][,4] <- M[[i-1]][,3]
        }
      }
      params <- matrix(unlist(lapply(M, get_params)),ncol=2,byrow=TRUE)
      resids <- lapply(M, get_resids)
      for (i in 1:frq){
        M[[i]][,3] <- resids[[i]]
      }
      if (sum(params - paramsx > 0.01) == 0 | n > 50) break
      paramsx <- params
    }
    # GENERATE WHITE NOISE FROM THE FITTED PARMA MODEL
    params <- cbind(params, unlist(lapply(resids,sd)))
    lenX <- years * frq * reps
    Xtr_ds <- vector("numeric", lenX)
    period <- rep(1:frq, lenX/frq)
    error <- rnorm(1, 0, params[1, 3])
    Xtr_ds[1] <- params[1, 1] *  tail(Qtr_ds,1) + params[1, 2] * tail(resids[[frq]],1) + error
    for (t in 2:length(Xtr_ds)){
      preverror <- error
      error <- rnorm(1, 0, params[period[t], 3])
      Xtr_ds[t] <- params[period[t], 1] * Xtr_ds[t - 1] + params[period[t], 2] * preverror + error
    }
    ## ADD SEASONS BACK
    Xtr <- as.vector(Xtr_ds * rep(tapply(Qtr, cyc, sd), lenX / frq)
                     + rep(tapply(Qtr, cyc, mean), lenX / frq))
    ## REVERSE THE LOG TRANSFORM
    a <- tapply(Q, cyc, getlogtransparam)
    mins <- tapply(Q, cyc, min)
    a[which(a>=mins)] <- mins[which(a>=mins)] - 1
    a <- rep(a, length(Xtr) / frq)
    X <- ts(matrix((exp(Xtr)) + a, ncol = reps), frequency = frq)
    if (missing(adjust)) adjust <- TRUE
    if (adjust==TRUE){
      X[which(X < 0)] <- 0
      X[which(X > 1.2 * max(Q))] <- 1.2 * max(Q)
    }
    params <- cbind(params, a[1:frq])
    colnames(params) <- c("phi", "theta", "st_dev", "log_trans_a")
  }
  if (missing(parameters) || parameters == FALSE){
    return(X)
  }else if (parameters == TRUE){
    results <- list(X, params)
    names(results) <- c("replicates", "parameters")
    return(results)
  }
}