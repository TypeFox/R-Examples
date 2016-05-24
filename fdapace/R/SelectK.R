#' Selects number of functional principal components for
#' given FPCA output and selection criteria
#'
#' @param fpcaObj A list containing FPCA related subjects returned by MakeFPCAResults().
#' @param criterion A string or positive integer specifying selection criterion for number of functional principal components, available options: 'FVE', 'AIC', 'BIC', or the specified number of components - default: 'AIC'
#' @param FVEthreshold A threshold percentage specified by user when using "FVE" as selection criterion: (0,1] - default: NULL
#' @param y A list of \emph{n} vectors containing the observed values for each individual - default: NULL
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y - default: NULL
#'
#' @return A list including the following two fields:
#' \item{k}{An integer indicating the selected number of components based on given criterion.}
#' \item{criterion}{The calculated criterion value for the selected number of components, i.e. FVE, AIC or BIC value, NULL for fixedK criterion.}
#'
#' @export

SelectK = function(fpcaObj, criterion = 'AIC', FVEthreshold = NULL, y = NULL, t = NULL){
  if(class(fpcaObj) != 'FPCA'){
    stop('Invalid Input: not a FPCA object!')
  }
  if(is.null(criterion)){
    stop('Invalid selection criterion. Selection criterion must not be NULL!')
  }
  if(!(criterion %in% c('FVE', 'AIC', 'BIC'))){
    if(is.numeric(criterion)){
      if(as.integer(criterion) != criterion || criterion <= 0){
        stop('Invalid selection criterion. To select fixed number of component, criterion needs to be a positive integer.')
      }
    } else {
      stop('Invalid selection criterion. Need to be one of "FVE", "AIC", "BIC" or a positive integer!')
    }
  }
  
  if(criterion %in% c('AIC','BIC')) {
    if(fpcaObj$optns$lean == TRUE && (is.null(y) || is.null(t))){
    stop("Option lean is TRUE, need input data y and measurement time list t to calculate log-likelihood.")
    }
    if(fpcaObj$optns$lean == FALSE){
      y <- fpcaObj$inputData$y
      t <- fpcaObj$inputData$t
    }
    if(criterion == 'AIC'){C = 2}
    else {C = log(length(y))}
    # FVE is not the selection criterion
    IC = rep(Inf, length(fpcaObj$lambda))
    for(i in 1:length(fpcaObj$lambda)){
      logliktemp = GetLogLik(fpcaObj, i, y = y, t = t)
      if(is.null(logliktemp)){
        stop('The covariance matrix of the estimated function is nearly singular! AIC or BIC is not applicable.')
      }
      IC[i] = logliktemp + C * i
      if(i > 1 && IC[i] > IC[i-1]){
        # cease whenever AIC/BIC stops decreasing
        return(list(k = i-1, criterion = IC[i-1]))
      }
      if(i == length(fpcaObj$lambda)){
        return(list(k = i, criterion = IC[i]))
      }
    }
    #if(criterion != 'FVE'){
    #  return(k = length(fpcaObj$lambda))
    #}
  } else if(criterion == 'FVE'){
    # select FVE based on cumFVE in fpcaObj and specified FVEthreshold
    if(is.null(FVEthreshold)){stop('Need to specify FVEthreshold to choose number of components via FVE.')}
    cumFVE = fpcaObj$cumFVE
    buff <- .Machine[['double.eps']] * 100
    return( list(k = min( which(cumFVE > FVEthreshold * 100 - buff) ), criterion = cumFVE[min(which(cumFVE > FVEthreshold * 100 - buff))]))
  } else { # fixed K is specified.
    if(criterion > length(fpcaObj$lambda)){
      stop("Specified number of components is more than available components.")
    }
    return(list(k = criterion, criterion = NULL))
  }
}
