#' Generic for fitting a Stochastic Mortality Model
#'
#' \code{fit} is a generic function for fitting Stochastic Mortality Models. 
#' The function invokes particular methods which depend on the class of the 
#' first argument. 
#'
#' \code{fit} is a generic function which means that new fitting strategies 
#' can be added for particular stochastic mortality models. See for instance 
#' \code{\link{fit.StMoMo}}.
#'
#' @param object an object used to select a method. Typically of class 
#' \code{StMoMo} or an extension of this class.
#'
#' @param ... arguments to be passed to or from other methods.
#'
#' @export
fit =  function(object, ...)
  UseMethod("fit")


#' Fit a Stochastic Mortality Model
#' 
#' Fit a Stochastic Mortality Model to a given data set. The fitting is done 
#' using package \code{gnm}.
#' 
#' Fitting is done using function \code{\link[gnm]{gnm}} within package 
#' \code{gnm}. This is equivalent to minimising (maximising) the deviance 
#' (log-likelihood) of the  model. Ages and years in the data should be of 
#' type numeric. Data points with zero exposure are assigned a zero weight 
#' and are ignored in the fitting process. Similarly, \code{NA} are assigned a
#' zero weight and ignored in the fitting process. Parameter estimates can be 
#' plotted using function \code{\link{plot.fitStMoMo}}.
#' 
#' @param object an object of class \code{"StMoMo"} defining the stochastic
#' mortality model.
#' @param Dxt matrix of deaths data.
#' @param Ext matrix of observed exposures of the same dimension of \code{Dxt}. 
#' @param ages vector of ages corresponding to rows of \code{Dxt} and \code{Ext}. 
#' @param years vector of years corresponding to rows of \code{Dxt} and \code{Ext}. 
#' @param ages.fit optional vector of ages to include in the fit. Must be a 
#' subset of \code{ages}. 
#' @param years.fit optional vector of years to include in the fit. Must be a 
#' subset of \code{years}. 
#' @param oxt optional matrix/vector or scalar of known offset to be used in fitting
#' the model. This can be used to specify any a priori known component to be added to 
#' the predictor during fitting. 
#' @param wxt optional matrix of 0-1 weights to be used in the fitting process. 
#' This can be used, for instance, to zero weight some cohorts in the data.
#' See \code{\link{genWeightMat}} which is a helper function for defining 
#' weighting matrices.
#' @param start.ax optional vector with starting values for \eqn{\alpha_x}.
#' @param start.bx optional matrix with starting values for \eqn{\beta_x^{(i)}}.
#' @param start.kt optional matrix with starting values for \eqn{\kappa_t^{(i)}}.
#' @param start.b0x optional vector with starting values for \eqn{\beta_x^{(0)}}.
#' @param start.gc optional vector with starting values for \eqn{\gamma_c}.
#' @param verbose a logical value. If \code{TRUE} progress indicators are 
#' printed as the model is fitted. Set \code{verbose = FALSE} to silent the 
#' fitting and avoid progress messages.
#' @param ... arguments to be passed to or from other methods. This can be 
#' used to control the fitting parameters of \code{gnm}. See 
#' \code{\link[gnm]{gnm}}.
#' 
#' @return A list with class \code{"fitStMoMo"} with components:
#'   
#'   \item{model}{ the object of class \code{"StMoMo"} defining the fitted 
#'   stochastic mortality model.}
#'   
#'   \item{ax}{ vector with the fitted values of the static age function 
#'   \eqn{\alpha_x}. If the model does not have a static age function or 
#'   failed to fit this is set to \code{NULL}.}
#'     
#'   \item{bx}{ matrix with the values of the period age-modulating functions 
#'   \eqn{\beta_x^{(i)}, i=1, ..., N}. If the \eqn{i}-th age-modulating 
#'   function is non-parametric (e.g. as in the Lee-Carter model) 
#'   \code{bx[, i]} contains the estimated values. If the model does not have 
#'   any age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{kt}{ matrix with the values of the fitted period indexes 
#'   \eqn{\kappa_t^{(i)}, i=1, ..., N}. \code{kt[i, ]} contains the estimated 
#'   values of the \eqn{i}-th period index. If the model does not have any 
#'   age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   \item{b0x}{ vector with the values of the cohort age-modulating function 
#'   \eqn{\beta_x^{(0)}}. If the age-modulating function is non-parametric 
#'   \code{b0x} contains the estimated values. If the model does not have a 
#'   cohort effect or failed to fit this is set to \code{NULL}.}
#'     
#'   \item{gc}{ vector with the fitted cohort index \eqn{\gamma_{c}}.
#'   If the model does not have a cohort effect or failed to fit this is set
#'   to \code{NULL}.}
#'   
#'   \item{Dxt}{ matrix of deaths used in the fitting.}
#'   
#'   \item{Ext}{ matrix of exposures used in the fitting.}
#'   
#'   \item{oxt}{ matrix of known offset values used in the fitting.}
#'   
#'   \item{wxt}{ matrix of 0-1 weights used in the fitting.}
#'   
#'   \item{ages}{ vector of ages used in the fitting.}
#'   
#'   \item{years}{ vector of years used in the fitting.}
#'   
#'   \item{cohorts}{ vector of cohorts used in the fitting.}
#'   
#'   \item{fittingModel}{ output from the call to \code{gnm} used to fit the 
#'   model. If the fitting failed to converge this is set to \code{NULL}.}
#'   
#'   \item{loglik}{ log-likelihood of the model. If the fitting failed to 
#'   converge this is set to \code{NULL}.}
#'   
#'   \item{deviance}{ deviance of the model. If the fitting failed to 
#'   converge this is set to \code{NULL}.}
#'  
#'   \item{npar}{ effective number of parameters in the model. If the fitting
#'   failed to converge this is set to \code{NULL}.}
#'    
#'    \item{nobs}{ number of observations in the model fit. If the fitting
#'    failed to converge this is set to \code{NULL}.}
#'
#'    \item{fail}{ \code{TRUE} if a model could not be fitted and 
#'    \code{FALSE} otherwise.}    
#'            
#'    \item{conv}{ \code{TRUE} if the model fitting converged and 
#'    \code{FALSE} if it didn't.}
#'     
#'  @seealso \code{\link{genWeightMat}}, \code{\link{plot.fitStMoMo}}
#'     
#' @examples    
#' 
#' # CBD model only to older ages
#' CBDfit <- fit(cbd(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               ages.fit = 55:89)
#' plot(CBDfit, parametricbx = FALSE)
#' 
#' # APC model weigthing out the 3 first and last cohorts
#' wxt <- genWeightMat(EWMaleData$ages,  EWMaleData$years, clip = 3)
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               wxt = wxt)
#' plot(APCfit, parametricbx = FALSE, nCol = 3)
#' 
#' # Set verbose = FALSE for silent fitting
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               wxt = wxt, verbose = FALSE)
#' \dontrun{
#' # Poisson Lee-Carter model with the static age function set to  
#' # the mean over time of the log-death rates
#' constLCfix_ax <- function(ax, bx, kt, b0x, gc, wxt, ages){  
#'   c1 <- sum(bx, na.rm = TRUE)
#'   bx <- bx / c1
#'   kt <- kt * c1  
#'   list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)  
#' }  
#' LCfix_ax <- StMoMo(link = "log", staticAgeFun = FALSE, 
#'                    periodAgeFun = "NP", constFun =  constLCfix_ax)
#' LCfix_axfit <- fit(LCfix_ax, Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'                    ages = EWMaleData$ages, years = EWMaleData$years, 
#'                    oxt = rowMeans(log(EWMaleData$Dxt / EWMaleData$Ext)))
#' plot(LCfix_axfit)
#' }
#' @export 
fit.StMoMo <- function(object, Dxt, Ext, ages = 1:nrow(Dxt), 
                       years = 1:ncol(Dxt), ages.fit = ages, 
                       years.fit = years, oxt = NULL, wxt = NULL, 
                       start.ax = NULL, start.bx = NULL, start.kt = NULL,
                       start.b0x = NULL, start.gc = NULL, verbose = TRUE, 
                       ...) {
  #Hack to remove notes in CRAN check
  x <- NULL
  w <- NULL
  
  # Construct fitting data
  
  nAges <- length(ages)
  nYears <- length(years)    
  Dxt <- as.matrix(Dxt)
  if (nrow(Dxt) != nAges ||  ncol(Dxt) != nYears) {
    stop( "Mismatch between the dimension of Dxt and the 
            number of years or ages")
  }
  rownames(Dxt) <- ages
  colnames(Dxt) <- years
  
  Ext <- as.matrix(Ext)
  if (nrow(Ext) != nAges ||  ncol(Ext) != nYears) {
    stop( "Mismatch between the dimension of Ext and the 
            number of years or ages")
  }
  rownames(Ext) <- ages
  colnames(Ext) <- years  
  
  #Extract the specific ages and years for fitting 
  if (length(ages.fit) != length(which(ages.fit %in% ages))) {
    stop( "ages.fit is not a subset of ages")
  }
  if (length(years.fit) != length(which(years.fit %in% years))) {
    stop( "years.fit is not a subset of years")
  }  
  Dxt <- Dxt[which(ages %in% ages.fit), which(years %in% years.fit)]
  Ext <- Ext[which(ages %in% ages.fit), which(years %in% years.fit)]  
  ages <- ages.fit
  years <- years.fit
  
  # Construct fitting data
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)  
  
  fitDataD<- (reshape2::melt(Dxt, value.name = "D", varnames = c("x", "t")))
  fitDataE<- (reshape2::melt(Ext, value.name = "E", varnames = c("x", "t")))
  fitData <- merge(fitDataD, fitDataE)      
  fitData <- transform(fitData, c = t - x)      
    
  if (is.null(oxt)) {
    oxt <- matrix(0, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years
    fitData$o <- 0      
  } else {
    oxt <- matrix(oxt, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years
    fitDataO <- (reshape2::melt(oxt, value.name = "o", 
                                varnames = c("x", "t")))
    fitData <- merge(fitData, fitDataO)       
  }   
  
  if (is.null(wxt)) {
    wxt <- matrix(1, nrow = nAges, ncol = nYears)
  } else {
    wxt <- as.matrix(wxt)
    if ( nrow(wxt) != nAges ||  ncol(wxt) != nYears) {
      stop( "Mismatch between the dimension of the weigth matrix wxt and the 
            number of fitting years or ages")
    }    
  } 
  rownames(wxt) <- ages
  colnames(wxt) <- years
  if (any(Ext <= 0)) { #Non-positive exposures
    indExt <- (Ext <= 0)
    wxt[indExt] <- 0
    warning(paste("StMoMo: ", sum(indExt), " data points have 
                  non-positive exposures and have been zero weighted\n", 
                  sep = ""))
  }
  if (any(is.na(Dxt / Ext))) { #Missing values
    indqxt <- is.na(Dxt / Ext)
    wxt[indqxt] <- 0
    warning(paste("StMoMo: ", sum(indqxt), 
                  " missing values which have been zero weighted\n", sep = ""))
  }
  
  fitDataW <- (reshape2::melt(wxt, value.name = "w", varnames = c("x", "t")))
  fitData <- merge(fitData, fitDataW)
  
  
  #data for age-parametric terms
  N <- object$N
  if (N > 0) {
    for (i in 1:N) {      
      if (is.function(object$periodAgeFun[[i]])) {
        f <- function(x) object$periodAgeFun[[i]](x, ages)
        fitData$B__ <- apply(as.array(fitData$x), MARGIN = 1, FUN = f)
        names(fitData)[names(fitData) == "B__"] <- paste("B", i, sep = "")
      }      
    }
  }
  
  #data for age-cohort term  
  if (is.function(object$cohortAgeFun)) {
    f <- function(x) object$cohortAgeFun(x, ages)
    fitData$B0 <- apply(as.array(fitData$x), MARGIN = 1, FUN = f)    
  }      
  
  
  # Remove from the data ages, years or cohorts with 0 weight 
  # This needs to be done as gnm can not handle this for non-parametric terms
  wxTemp <- aggregate(data = fitData, w ~ x, FUN = sum)  
  zeroWeigthAges <- as.character(wxTemp$x[which((wxTemp$w <= 0))])
  wtTemp <- aggregate(data = fitData, w ~ t, FUN = sum)  
  zeroWeigthYears <- as.character(wtTemp$t[which((wtTemp$w <= 0))])  
  wcTemp <- aggregate(data = fitData, w ~ c, FUN = sum)  
  zeroWeigthCohorts <- as.character(wcTemp$c[which((wcTemp$w <= 0))])      
  fitData <- subset(fitData, w > 0)
  
  if (verbose) {
    if (length(zeroWeigthAges) > 0) {
      cat("StMoMo: The following ages have been zero weigthed:", 
          zeroWeigthAges, "\n")
    }
    if (length(zeroWeigthYears) > 0) {
      cat("StMoMo: The following years have been zero weigthed:", 
          zeroWeigthYears, "\n")
    }
    if (length(zeroWeigthCohorts) > 0) {
      cat("StMoMo: The following cohorts have been zero weigthed:", 
          zeroWeigthCohorts, "\n")
    }    
  }  
  
  ### Set starting values
  if (object$link == "logit") {        
    coefNames <- gnm(formula = as.formula(object$gnmFormula), data = fitData,
                     family = binomial(link = "logit"), 
                     weights = fitData$E * fitData$w, start = start, 
                     method = "coefNames")
  } else {        
    coefNames <- gnm(formula = as.formula(object$gnmFormula), data = fitData,
                     family = poisson(link = "log"), 
                     weights = fitData$w, start = start, 
                     method = "coefNames")
  }  
  if (is.null(start.ax)) start.ax <- rep(NA, nAges)
  if (is.null(start.bx)) start.bx <- array(NA, c(nAges, N))
  if (is.null(start.kt)) start.kt <- array(NA, c(N, nYears))
  if (is.null(start.b0x)) start.b0x <- rep(NA, nAges)
  if (is.null(start.gc)) start.gc <- rep(NA, nCohorts)
  startCoef <- processStartValues(object, coefNames, start.ax, start.bx, 
                                  start.kt, start.b0x, start.gc, ages, 
                                  years, cohorts)
  
  
  ### Fit using gnm  
  if (verbose) cat("StMoMo: Start fitting with gnm\n")  
  if (object$link == "logit") {        
    fittingModel <- gnm(formula = as.formula(object$gnmFormula), 
                        data = fitData, family = binomial(link = "logit"),
                        weights = fitData$E * fitData$w, start = startCoef, 
                        verbose = verbose, ...)
  } else {        
    fittingModel <- gnm(formula = as.formula(object$gnmFormula), 
                        data = fitData, family= poisson(link = "log"), 
                        weights = fitData$w, start = startCoef, 
                        verbose = verbose, ...)
  }  
  if (verbose) 
    cat("StMoMo: Finish fitting with gnm\n")  
  fail <- is.null(fittingModel$conv)
  if (fail) {
    conv <- FALSE
    warning("StMoMo: The model fitting failed and no model could be estimated.\n")      
    out <- list(model = object, ax = NULL, bx = NULL, kt = NULL, b0x = NULL, 
                gc = NULL, Dxt = Dxt, Ext = Ext, oxt = oxt , wxt = wxt, 
                ages = ages, years = years, cohorts = cohorts,
                fittingModel = fittingModel, loglik = NULL,
                deviance = NULL, npar = NULL, nobs = NULL, conv = conv, 
                fail = fail, fitData = fitData)          
    class(out) <- "fitStMoMo"
    return(out)    
  }
  else {
    conv <- fittingModel$conv[1]  
    if (!conv) warning("StMoMo: The model fitting didn't converge.\n")  
  }  
  
  ### Extract coefficients    
  fittedCoef <- extractCoefficientsFromGnm(object, coef(fittingModel), ages, 
                                           years, cohorts, zeroWeigthAges, 
                                           zeroWeigthYears, zeroWeigthCohorts)  
  
  ### Apply identifiability constraints
  oldLinkPred <- predictLink(fittedCoef$ax, fittedCoef$bx, fittedCoef$kt, 
                             fittedCoef$b0x, fittedCoef$gc, oxt, ages, years)
  constPar<- object$constFun(fittedCoef$ax, fittedCoef$bx, fittedCoef$kt, 
                             fittedCoef$b0x, fittedCoef$gc, wxt, ages)
  ax <- constPar$ax
  bx <- constPar$bx
  kt <- constPar$kt
  b0x <- constPar$b0x
  gc <- constPar$gc
  
  # Check if NA where introduced by th indentifiability constraints
  NAinNewParameters <- (!is.null(fittedCoef$ax) && 
                          (is.null(ax) || anyNA(ax[!is.na(fittedCoef$ax)]))) || 
                       (!is.null(fittedCoef$bx) && 
                          (is.null(bx) || anyNA(bx[!is.na(fittedCoef$bx)]))) ||
                       (!is.null(fittedCoef$kt) && 
                          (is.null(kt) || anyNA(kt[!is.na(fittedCoef$kt)]))) ||     
                       (!is.null(fittedCoef$b0x) && 
                          (is.null(b0x) || anyNA(b0x[!is.na(fittedCoef$b0x)]))) ||   
                       (!is.null(fittedCoef$gc) && 
                          (is.null(gc) || anyNA(gc[!is.na(fittedCoef$gc)])))
                       
  if (NAinNewParameters) {
    stop("The parameter trasnformation function transformed some parameters into NA or NULL.
         Check the 'constFun' argument of StMoMo.\n")
  }
  
  # Check if the transformation preserves the link.
  newLinkPred <- predictLink(ax, bx, kt, b0x, gc, oxt, ages, years)    
  trasError <- max(abs((newLinkPred[wxt != 0] - oldLinkPred[wxt !=  0]) / oldLinkPred[wxt !=  0]), 
                   na.rm = TRUE)
  if (trasError > 1e-10){
    stop("The parameter transformation function does not preserve the fitted rates.
          Check the 'constFun' argument of StMoMo.\n")
  }
  
  
  ### Preparet the output  
  
  # Compute log-like like-lihood and the deviance  
  if (object$link == "logit") {
    temp <- exp(newLinkPred)
    qhatxt <- temp / (temp + 1)
    loglik <- computeLogLikBinomial(obs = Dxt / Ext, fit = qhatxt, 
                                    exposure = Ext, weight = wxt)
    deviance <- computeDevianceBinomial(obs = Dxt / Ext, fit = qhatxt,
                                        exposure = Ext, weight = wxt)      
  } else if (object$link == "log") {
    Dhatxt <- exp(newLinkPred) * Ext
    loglik <- computeLogLikPoisson(obs = Dxt, fit = Dhatxt, weight = wxt)
    deviance <- computeDeviancePoisson(obs = Dxt, fit = Dhatxt, weight = wxt)
  }
  
  out <- list(model = object, ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc, 
              Dxt = Dxt, Ext = Ext, oxt = oxt , wxt = wxt, ages = ages, 
              years = years, cohorts = cohorts, fittingModel = fittingModel, 
              loglik = loglik, deviance = deviance,   
              npar = fittingModel$rank[1], nobs = nobs(fittingModel), 
              conv = conv, fail = fail, call  = match.call())          
  class(out) <- "fitStMoMo"
  out    
}


#' Extract the model coefficient
#' 
#' Extract the model coefficient of a stochastic mortality model from
#' a gnm fit of the model
#' 
#' @param object an object of class \code{"StMoMo"} defining the stochastic
#' mortality model.
#' @param coefGnmModel fitted coefficient from a gnm model fit.
#' @param ages ages in the fitting data.
#' @param years years in the fitting data.
#' @param cohorts cohorts in the fitting data.
#' @param zeroWeigthAges character vector of years whose parameters cannot 
#' be estimated because all data is zero weighted
#' @param zeroWeigthYears character vector of years whose parameters cannot 
#' be estimated because all data is zero weighted
#' @param zeroWeigthCohorts character vector of cohort whose parameters 
#' cannot be estimated because all data is zero weighted
#' 
#' @details Weigth vectors wx, wx, wc are used to identify parameters that
#' cannot be estimated because all the data is weighted out.
#' 
#' @return A list with the model parameters, ax, bx, kt, b0x, gc
#' @keywords internal
extractCoefficientsFromGnm <- function(object,coefGnmModel, ages, years, 
                                       cohorts, zeroWeigthAges, 
                                       zeroWeigthYears, zeroWeigthCohorts) {
  
  nAges <- length(ages)
  nYears <- length(years)  
  nCohorts <- length(cohorts)
  N <- object$N
  
  # age term
  if (object$staticAgeFun == TRUE) {
    axTemp <- coefGnmModel[grep(pattern = "^factor[(]x[)]", 
                                names(coefGnmModel))]    
    names(axTemp) <- sub(pattern = "factor[(]x[)]", replacement = "" , 
                         x = names(axTemp))
    ax <- rep(NA,nAges)
    names(ax) <- ages
    ax[names(axTemp)] <- axTemp
    ax[zeroWeigthAges] <- NA
  }
  else{
    ax <- NULL
  }
  
  # age-period terms
  bx <- NULL
  kt <- NULL
  if (N > 0) {    
    bx <- array(1, dim = c(nAges, N), dimnames = list(ages, 1:N))
    kt <- array(0, dim = c(N, nYears), dimnames = list(1:N, years))
    for (i in 1:N) {
      if (is.function(object$periodAgeFun[[i]])) {
        #bx
        f <- function(x) object$periodAgeFun[[i]](x, ages)
        bx[, i] <- apply(as.array(ages), MARGIN = 1, FUN = f)
        
        #kt
        pattern1 <- paste("^factor[(]t[)][-]?[[:digit:]]+:B", i, sep = "")
        pattern2 <- paste("^B",i,":factor[(]t[)][-]?[[:digit:]]", sep = "")
        ktTemp <- coefGnmModel[c(grep(pattern = pattern1, names(coefGnmModel)),
                                 grep(pattern = pattern2, names(coefGnmModel)))] 
        names(ktTemp) <- sub(pattern = "factor[(]t[)]", replacement = "" ,
                             x = names(ktTemp))
        names(ktTemp) <- sub(pattern = paste("B", i, sep = ""), 
                             replacement = "" , x = names(ktTemp))
        names(ktTemp) <- sub(pattern = ":", replacement = "" , 
                             x = names(ktTemp))
        ktTemp[is.na(ktTemp)] <- 0
        kt[i, names(ktTemp)] <- ktTemp 
        kt[i, zeroWeigthYears] <- NA 
      } else if (object$periodAgeFun[[i]] == "1") {        
        #kt
        pattern <- "^factor[(]t[)][-]?[[:digit:]]+$"
        ktTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))] 
        names(ktTemp) <- sub(pattern = "factor[(]t[)]", replacement = "" , 
                             x = names(ktTemp))
        ktTemp[is.na(ktTemp)] <- 0
        kt[i, names(ktTemp)] <- ktTemp 
        kt[i, zeroWeigthYears] <- NA 
      } else {  # Non-Parametric terms
        ind <- i:N
        inst <- 1 + sum(object$periodAgeFun[-ind] == object$periodAgeFun[[i]])         
        pattern <- paste("Mult[(]., factor[(]t[)], inst = ",inst,
                         "[)].factor[(]x[)]", sep = "")
        bxTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))]
        names(bxTemp) <- sub(pattern = pattern, replacement = "" , 
                             x = names(bxTemp))
        bx[names(bxTemp), i] <- bxTemp
        bx[zeroWeigthAges, i] <- NA 
        #kt
        pattern <- paste("Mult[(]factor[(]x[)], ., inst = ", inst, 
                         "[)].factor[(]t[)]", sep = "")
        ktTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))] 
        names(ktTemp) <- sub(pattern = pattern, replacement = "" , 
                             x = names(ktTemp))
        ktTemp[is.na(ktTemp)] <- 0
        kt[i, names(ktTemp)] <- ktTemp      
        kt[i, zeroWeigthYears] <- NA 
      }      
      
    }
    
  }
  
  #age-cohort term
  if ( !is.null(object$cohortAgeFun) == TRUE ) {
    b0x <- rep(1, nAges)
    names(b0x) <- ages    
    gc <- rep(0, nCohorts)
    names(gc) <- cohorts
    if (is.function(object$cohortAgeFun)) {
      #b0x
      f <- function(x) object$cohortAgeFun(x, ages)
      b0x <- apply(as.array(ages), MARGIN = 1, FUN = f)
      names(b0x) <- ages        
      #gc
      pattern1 <- "^factor[(]c[)][-]?[[:digit:]]+:B0"
      pattern2 <- "^B0:factor[(]c[)][-]?[[:digit:]]"
      gcTemp <- coefGnmModel[c(grep(pattern = pattern1, names(coefGnmModel)),
                               grep(pattern = pattern2, names(coefGnmModel)))] 
      names(gcTemp) <- sub(pattern = "factor[(]c[)]", replacement = "" , 
                           x = names(gcTemp))
      names(gcTemp) <- sub(pattern = "B0", replacement = "" , x = names(gcTemp))
      names(gcTemp) <- sub(pattern = ":", replacement = "" , x = names(gcTemp))
      gcTemp[is.na(gcTemp)] <- 0      
      gc[names(gcTemp)] <- gcTemp
      gc[zeroWeigthCohorts] <- NA
      
    } else if (object$cohortAgeFun == "1") {        
      #gc
      pattern <- "^factor[(]c[)][-]?[[:digit:]]+$"
      gcTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))] 
      names(gcTemp) <- sub(pattern = "factor[(]c[)]", replacement = "" , 
                           x = names(gcTemp))
      gcTemp[is.na(gcTemp)] <- 0
      gc[names(gcTemp)]<-gcTemp 
      gc[zeroWeigthCohorts] <- NA
      
    } else {  # Non-Parametric terms
      #b0x
      pattern <- "Mult[(]., factor[(]c[)][)].factor[(]x[)]"
      b0xTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))]
      names(b0xTemp) <- sub(pattern = pattern, replacement = "" , 
                            x = names(b0xTemp))
      b0x[names(b0xTemp)] <- b0xTemp
      b0x[zeroWeigthAges] <- NA
      #gc
      pattern <- "Mult[(]factor[(]x[)], .[)].factor[(]c[)]"
      gcTemp <- coefGnmModel[grep(pattern = pattern, names(coefGnmModel))] 
      names(gcTemp) <- sub(pattern = pattern, replacement = "" , 
                           x = names(gcTemp))
      gcTemp[is.na(gcTemp)] <- 0
      gc[names(gcTemp)] <- gcTemp
      gc[zeroWeigthCohorts] <- NA
      
    }
  } else {
    b0x <- NULL
    gc <- NULL
  }   
  
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}

#' Process the initial parameter supplied by the user
#' 
#' Convert the initial parameters supplied by the user into parameters 
#' suitable for use within \code{gnm}
#' 
#' @param object an object of class \code{"StMoMo"} defining the stochastic
#' mortality model.
#' @param coefNames name of the coefficients in the gnm model
#' @param ax vector with starting values for \eqn{\alpha[x]}
#' @param bx matrix with starting values for \eqn{\beta[x]}
#' @param kt matrix with starting values for \eqn{\kappa[t]}
#' @param b0x vector with starting values for \eqn{\beta_x^{(0)}}
#' @param gc vector with starting values for \eqn{\gamma[c]}
#' @param ages ages in the fitting data.
#' @param years years in the fitting data.
#' @param cohorts cohorts in the fitting data.
#' 
#' @return a vector of intitial parameter estimates
#' 
#' @keywords internal
processStartValues <- function(object, coefNames, ax, bx, kt, b0x, gc, 
                               ages, years, cohorts) {
  
  startCoef <- rep(NA, length(coefNames))
  names(startCoef) <- coefNames
  nAges <- length(ages)
  nYears <- length(years)  
  nCohorts <- length(cohorts)
  N <- object$N
  
  # age term
  if (object$staticAgeFun == TRUE) { 
    if (length(ax) != nAges) {
      stop( "Mismatch between the number of ages and start.ax")
    }
    names(ax) <- ages
    for (x in as.character(ages)) {
      ind <- which(coefNames == paste("factor(x)", x, sep = ""))      
      if (length(ind) > 0) startCoef[ind] <- ax[x]      
    }   
  }
  
  # age-period terms  
  if (N > 0) {    
    if (nrow(bx) != nAges) {
      stop( "Mismatch between the number of ages and start.bx.")
    }
    if (ncol(bx) != N) {
      stop( "Mismatch between the number of age/period terms and start.bx.")
    }
    if (ncol(kt) != nYears) {
      stop( "Mismatch between the number of years and start.kt.")
    }
    if (nrow(kt) != N) {
      stop( "Mismatch between the number of age/period terms and start.kt.")
    }
    
    dimnames(bx) <- list(ages, 1:N)    
    dimnames(kt) <- list(1:N, years)
    for (i in 1:N) {
      if (is.function(object$periodAgeFun[[i]])) {        
        #kt
        for (t in as.character(years)) {          
          ind1 <- which(coefNames == paste("factor(t)", t,":B", i, sep = ""))
          ind2 <- which(coefNames == paste("B", i,":factor(t)", t, sep = ""))
          ind <- c(ind1, ind2)
          if (length(ind) > 0) startCoef[ind] <- kt[i, t]        
        }
      } else if (object$periodAgeFun[[i]] == "1") {        
        #kt
        for (t in as.character(years)) {
          ind <- which(coefNames == paste("factor(t)", t, sep = ""))                
          if (length(ind) > 0) startCoef[ind] <- kt[i, t]                
        }
      } else {  # Non-Parametric terms
        
        indInst <- i:N
        inst <- 1 + sum(object$periodAgeFun[-indInst] == object$periodAgeFun[[i]])                
        #bx
        for (x in as.character(ages)) {          
          ind <- which(coefNames == paste("Mult(., factor(t), inst = ", inst,
                                          ").factor(x)", x, sep = ""))
          if (length(ind) > 0) startCoef[ind] <- bx[x, i]
        }        
        #kt
        for (t in as.character(years)) {
          ind <- which(coefNames == paste("Mult(factor(x), ., inst = ", inst, 
                                          ").factor(t)", t, sep = ""))
          if (length(ind) > 0) startCoef[ind] <- kt[i, t]                
        }        
      }            
    }
  }  
  
  #age-cohort term
  names(b0x) <- ages
  names(gc) <- cohorts
  if ( !is.null(object$cohortAgeFun) == TRUE ) {
    if (length(b0x) != nAges && object$cohortAgeFun == "NP") {
      stop( "Mismatch between the number of ages and start.b0x.")
    }
    if (length(gc) != nCohorts) {
      stop( "Mismatch between the number of cohorts and start.gc.")
    }    
    if (is.function(object$cohortAgeFun)) {
      #gc     
      for (k in as.character(cohorts)) {          
        ind1 <- which(coefNames == paste("factor(c)", k,":B0", sep = ""))
        ind2 <- which(coefNames == paste("B0:factor(c)", k, sep = ""))
        ind <- c(ind1, ind2)               
        if (length(ind) > 0) startCoef[ind] <- gc[k]        
      }      
    } else if (object$cohortAgeFun == "1") {        
      #gc     
      for (k in as.character(cohorts)) {          
        ind <- which(coefNames == paste("factor(c)", k, sep = ""))
        if (length(ind) > 0) startCoef[ind] <- gc[k]        
      }      
    } else {  # Non-Parametric terms
      #b0x
      for (x in as.character(ages)) {          
        ind <- which(coefNames == paste("Mult(., factor(c)).factor(x)", x, 
                                        sep = ""))
        if (length(ind) > 0) startCoef[ind] <- b0x[x]
      }        
      #gc
      for (k in as.character(cohorts)) {
        ind <- which(coefNames == paste("Mult(factor(x), .).factor(c)", k, 
                                        sep = ""))
        if (length(ind) > 0) startCoef[ind] <- gc[k]                
      }      
    }
  }
  startCoef
}


#' Print an object of class \code{"fitStMoMo"}
#' 
#' \code{print} method for class \code{"fitStMoMo"}. 
#' @usage 
#' \method{print}{fitStMoMo}(x, ...)
#' @param x an object of class \code{"fitStMoMo"}.
#' @param ... arguments to be passed to or from other methods.
#' @export 
#' @method print fitStMoMo
print.fitStMoMo <- function(x, ...) {
  cat("Stochastic Mortality Model fit")
  cat(paste("\nCall:", deparse(x$call)))
  cat("\n\n")
  print(x$model)  
  cat(paste("\n\nYears in fit:", min(x$years), "-", max(x$years)))
  cat(paste("\nAges in fit:", min(x$ages), "-", max(x$ages), "\n"))
  
  cat(paste("\nLog-likelihood: ", round(x$loglik[1], 2)))
  cat(paste("\nDeviance: ", round(x$deviance[1], 2)))
  cat(paste("\nNumber of parameters: ", x$npar))  
}

