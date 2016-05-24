#' Extended multiplicative signal correction (EMSC)
#'
#' Performs model-based background correction and normalisation of spectra. EMSC handles
#' variations in scaling, polynomial baselines and interferents. Parameters for corrections
#' are stored for further analysis, and spectra are corrected accordingly.
#'
#' This is the main EMSC function performing all calculations. It can be run with
#' no parameters (defaults are used), with a predefined EMSC model object or with
#' parameters that are passed on to the EMSC model building function \code{\link{EMSC_model}}.
#'
#' @param X \code{matrix} containing spectra as rows.
#' @param model an EMSC model to use instead of the other parameters.
#' @param ... named model parameters for EMSC_model.
#'
#' @return An object of class EMSC is returned. This contains:
#' \itemize{
#'  \item{\code{corrected}:}{ \code{matrix} of corrected spectra.}
#'  \item{\code{parameters}:}{ \code{matrix} of fitted parameter values.}
#'  \item{\code{model}:}{ object containing input all input parameters.}
#'  \item{\code{X}:}{ original data.}
#' }
#'
#' @seealso \code{\link{EMSC_model}} \code{\link{predict.EMSC}} \code{\link{plot.EMSC}}
#' 
#' @references H. Martens, E. Stark, Extended multiplicative signal correction and spectral
#'  interference subtraction: new preprocessing methods for near infrared spectroscopy.
#'  J Pharm Biomed Anal. 1991; 9(8):625-35.
#'
#' @examples
#' data(milk)
#' Raman      <- milk$Raman[, 850:3300]
#' EMSC.basic <- EMSC(Raman)
#' EMSC.poly6 <- EMSC(Raman, degree = 6)
#' EMSC.rep   <- EMSC(Raman, degree = 6, reference = Raman[30, ],
#'                    replicates = milk$replicates)
#'
#' \dontrun{
#' old.par  <- par(mfrow = c(2,2), mar = c(4,4,1,1))
#' xlim     <- rev(as.numeric(range(colnames(Raman))))
#' matplot(colnames(Raman), t(Raman), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Raw spectra')
#' matplot(colnames(Raman), t(EMSC.basic$corrected), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Corrected (basic)')
#' matplot(colnames(Raman), t(EMSC.poly6$corrected), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Corrected (6th degree polynomial)')
#' matplot(colnames(Raman), t(EMSC.rep$corrected),   type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', 
#'         xlab = 'Corrected (reference = spec. #30, replicate correction (90%))')
#' par(old.par)
#' }
#'
#' @importFrom pracma mldivide
#' @export
EMSC <- function(X, model = NULL, ...){
  n <- dim(X)[1]

  # Strip X
  X <- as.matrix(X)
  X <- unclass(X)
  
  # Make EMSC model if not supplied
  if(is.null(model)){
    mf <- match.call(expand.dots = TRUE)
    if(is.null(x <- mf$x)){
      x <- X
    }
    model <- EMSC_model(x = X, ...)
  }

  # Extract from model object
  terms <- model$sizes[2] + 2 + model$sizes[3] + model$sizes[5] # Polynomial + interferents + replicate model

  # Apply weights (if included)
  if(!is.null(model$weights)){
    Xw     <- X * rep(model$weights, each = n)
    modelw <- model$model * rep(model$weights, each = n)
  } else {
    Xw     <- X
    modelw <- model$model
  }

  # Handling of 0 objects
  PP <- apply(Xw==0, 1, all)

  # Perform parameter estimation
  if(any(PP)){
    Parameters0  <- mldivide(t(modelw), t(Xw[!PP,]), pinv = FALSE)

    k <- 2:terms
    Corrected0 <- Xw[!PP,] - crossprod(Parameters0[k,, drop = FALSE], modelw[k,, drop = FALSE])

    corrected        <- Xw
    corrected[!PP,]  <- Corrected0
    parameters       <- matrix(0, nrow(Parameters0), nrow(Xw))
    parameters[,!PP] <- Parameters0

  } else {
    Parameters0  <- mldivide(t(modelw), t(Xw), pinv = FALSE)

    k <- 2:terms
    Corrected0 <- Xw - crossprod(Parameters0[k,, drop = FALSE], modelw[k,, drop = FALSE])
    Corrected0 <- Corrected0 * 1./Parameters0[1,] # correct multipl. eff.

    corrected  <- Corrected0
    parameters <- Parameters0
  }

  # Return
  object <- list(corrected = corrected, parameters = parameters, model = model, X = X)
  class(object) <- c("EMSC", "list")
  object
}


#' Model object for extended multiplicative signal correction (EMSC)
#'
#' Sets up an EMSC model to be applied to one or more set of spectra.
#'
#' @param x \code{numeric} vector containing abcissas of spectra to be corrected or matrix to be 
#' corrected with/without names colnames.
#' @param reference \code{numeric} vector containing the reference spectrum.
#' @param degree \code{integer} giving the polynomial degree of the baseline; 0 or higher, default = 2.
#' @param interferent \code{numeric} vector containing a spectral component to remove.
#' @param constituent \code{numeric} vector containing a spectral component to include.
#' @param weights \code{numeric} vector of abcissas weights.
#' @param replicates optional \code{vector} which identifies replicates. Default = NULL, 
#' meaning no replicate correctio will be performed.
#' @param rep_corr proportion of variance or number of subspace components in replicate space (default = 0.9).
#'
#' @return An EMSC model is returned containing all parameters.
#'
#' @seealso \code{\link{EMSC}} \code{\link{predict.EMSC}} \code{\link{plot.EMSC}}
#' @export
EMSC_model <- function(x, reference = NA, degree = 2,
                       interferent = NULL, constituent = NULL, weights = NULL,
                       replicates = NULL, rep_corr = 0.9){

  # Handle abcissas
  if(is.data.frame(x))
    x <- as.matrix(x)
  if(is.matrix(x)){
    if(!is.ordered(abcissas <- colnames(x))){
      abcissas <- 1:dim(x)[2]
    }
    if(length(reference) == 1 && is.na(reference))
      reference <- colMeans(x)
  } else {
    if(!is.ordered(abcissas <- names(x))){
      if(!is.ordered(abcissas <- x)){
        abcissas <- 1:length(x)
      }
    }
  }

  # Sizes
  n.i <- ifelse(is.matrix(interferent), nrow(interferent), ifelse(is.null(interferent),0, 1))
  n.c <- ifelse(is.matrix(constituent), nrow(constituent), ifelse(is.null(constituent),0, 1))
  p <- length(abcissas)

  Start   <- abcissas[1]
  End     <- abcissas[p]
  C  <- 0.5*(Start+End)
  M  <- 2.0/(Start-End)

  # Construct polynomials
  if(is.null(degree)){
    degree <- -1
  }
  model <- matrix(0, degree+2, p)
  model[1,] <- reference
  mod_names <- "Reference"
  if(!is.null(degree)){ # Contains at least baseline
    model[2,] <- 1
    mod_names <- c(mod_names, "Baseline")
    if(degree > 0){ # Contains polynomial
      for(i in 1:degree){
        model[i+2,] <- (M*(abcissas-C))^i
        mod_names <- c(mod_names, paste("Degree", i))
      }
      model[3,] <- M*(Start-abcissas)-1
    }
  }
  
  # Orthogonalize model spectra (except reference) to avoid overlap
  if(degree > 2){
    model[-1, ] <- t(qr.Q(qr(t(model[-1, , drop = FALSE])))[,1:(degree+1)])
    model[ 2, ] <- abs(model[2, ])
  }
  
  # Add interferents and constituents
  model <- rbind(model, interferent, constituent)
  if(!is.null(interferent)){
    mod_names <- c(mod_names, paste("Interferent", 1:n.i))
  }
  if(!is.null(constituent)){
    mod_names <- c(mod_names, paste("Constituent", 1:n.c))
  }

  # Add dimnames
  dimnames(model) <- list(term = mod_names, abcissas = abcissas)
  
  # Handle replicate correction
  if(!is.null(replicates) && length(replicates) > 1){
    M <- calculateReplicateModel(x, replicates, rep_corr,
                                     list(model = model, 
                                          abcissas = abcissas, 
                                          sizes = c(p, degree, n.i, n.c, 0), 
                                          weights = weights))
    return(M)
  } else { # ... or return without
    return(list(model = model, 
                abcissas = abcissas, 
                sizes = c(p, degree, n.i, n.c, 0), 
                weights = weights))
  }
}


# Not exported calculation of replicate vectors
calculateReplicateModel <- function(x, replicates, rep_corr, model){
  # Remove singletons
  n_reps <- table(replicates)
  if(any(n_reps == 1)){
    out <- numeric(0)
    for(i in unique(replicates)){
      if(sum(ind <- replicates == i) == 1){
        out <- c(out, which(ind))
      }
    }
    x_rep <- x[-out, , drop = FALSE]
    replicates <- replicates[-out]
  } else {
    x_rep <- x
  }
  
  # Save reference
  reference <- model$model[1,]
  
  # Loop over replicate sets
  uniqueReplicates <- unique(replicates)
  nreps    <- length(uniqueReplicates)
  meanCorr <- matrix(0, nreps, model$sizes[1])
  resCorr  <- matrix(0, 0, model$sizes[1])
  j <- 0
  for(i in uniqueReplicates){
    # Correct current replicate set
    reps <- x[replicates == i, ]
    j <- j+1
    model$model[1,] <- colMeans(x) 
    repsCorrected <- EMSC(reps, model = model)
    
    # Compute replicate means and corrected residuals
    meanCorr[j,] <- colMeans(repsCorrected$corrected)
    resCorr      <- rbind(resCorr, repsCorrected$corrected - rep(meanCorr[j,], each = dim(reps)[1]))
  }
  
  # Weight before PCA
  if(!is.null(model$weights)){
    resCorr <- resCorr * rep(model$weights, each = dim(resCorr)[1])
  }
  usv <- svd(resCorr)
  
  # Calculate explained variance
  ExpVar <- usv$d^2 / sum(usv$d^2)
  ExpVarCum <- cumsum(ExpVar)

  # Select number of components
  if(rep_corr < 1){
    rep_corr <- match(TRUE, ExpVarCum >= rep_corr)
  }
  
  # Insert replicate model as interferents
  rep_vecs <- usv$v[,1:rep_corr, drop = FALSE]
  colnames(rep_vecs) <- paste("Rep.vec.",1:rep_corr)
  if(model$sizes[4] == 0){
    model$model <- rbind(model$model, t(rep_vecs))
  } else {
    degInt <- model$sizes[2]+model$sizes[3]+2
    model$model <- rbind(model$model[1:degInt],
                         t(rep_vecs), 
                         model$model[(degInt+1):(degInt+model$sizes[3])])
  }
  model$model[1,] <- reference
  model$sizes[5]  <- rep_corr

  # Save replicate model
  model$rep <- list(PCA    = usv,
                    ExpVar = ExpVarCum,
                    ncomp  = rep_corr,
                    deviations = resCorr)
  model
}
