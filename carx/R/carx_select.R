#' Select a \code{carx} model by the AIC 
#'
#' This function selects the \code{carx} model which minimizes the AIC among a set of \code{carx} models
#'  defined by a set of  formulas or a list of regression formulas with a maximal AR order.
#' The model specification is supplied by \code{formulas} which can be either a formula or a list of formulas.
#' For each formula, the function will estimate the \code{carx} models with the AR order
#'  from 1 to \code{max.ar} inclusive.
#' If \code{detect.outlier=TRUE}, outlier detection will be performed for each combination of model 
#' formula and AR order.
#' The function returns a \code{list} which consists of: 1)  \code{aicMat} which is a matrix of AIC
#' values where 
#' each row contains the AICs of the model given by a specific regression formula with the AR order 
#' ranging from 1 to \code{mar.ar} 
#' (after incorporation of any found outlier if outlier detection if enabled), and 
#' 2) \code{fitted} which is the fitted object of the selected model.

#' @param formulas a regression formula or a  list of regression formulas.
#' @param data a \code{CenTS} object containing the data and censored information.
#' @param max.ar the maximal AR order.
#' @param detect.outlier logical to specify whether outlier detection is performed (and incorporating in 
#' the \code{carx} model any found additive outliers) before computing the AIC for a model. Default = \code{FALSE}.
#' @param ... other arguments to be supplied, if not null, it will be called with the selected model and data. 
#' Examples include \code{CI.compute=TRUE}, which will cause the function to re-estimate the selected model with the confidence 
#' intervals computed, as in the selection part, no confidence interval is  computed.
#' @export
#' @return a list consisting of:
#' \itemize{
#' \item{\code{fitted}}{ the fitted CARX object of the model with the smallest AIC.}
#' \item{\code{aicMat}}{ the matrix of AIC where rows correspond to the model formulas and columns correspond to the AR orders.}
#' }
#' @examples
#' dataSim <- carxSimCenTS(nObs=100)
#' fmls <- list(M1=y~X1,M2=y~X1+X2,M3=y~X1+X2-1)
#' \dontrun{cs = carxSelect(y~X1,max.ar=3,data=dataSim)}
#' \dontrun{cs = carxSelect(formulas=fmls,max.ar=3,data=dataSim)}
#' \dontrun{
#'   #To compute confidence intervals for the selected model, call with CI.compute=TRUE.
#'   cs = carxSelect(formulas=fmls,max.ar=3,data=dataSim,CI.compute=TRUE)
#' }
#'
carxSelect <- function(formulas, max.ar, data=list(), detect.outlier=F
                       #,verbose=FALSE
                       ,...)
{
  dotsArgs = list(...)
  if(typeof(formulas) == 'language')
    formulas <- list(M1=formulas)
  if(is.null(formulas))
    names(formulas) <- paste0("M",1:length(formulas))

  nModel <- length(formulas)
  aics <- matrix(nrow=nModel,ncol=max.ar)
  rownames(aics) <- names(formulas)
  colnames(aics) <- paste0("AR",seq(1,max.ar))
  saic0 <- NULL
  m0 <- NULL

  i  <- 0
  for(f in formulas)
  {
    i <- i + 1
    for( p in 1:max.ar)
    {
      #message(paste("Trying model fomula:",deparse(f), ", with AR order", p))
      toPass = c(list(formula=f,data=data,p=p),dotsArgs)
      toPass$CI.compute=FALSE #to save time
      #tmp <- carx(f,data=data, p=p, CI.compute=FALSE,...)
      tmp <- do.call(carx,toPass)
      if(detect.outlier)
        tmp <- outlier(tmp)
      a <- AIC(tmp)
      #if(verbose) message(paste0("Model formula:", deparse(formula(tmp)), ", AR order:",tmp$p, ", AIC: ", round(a,digits=4)))
      aics[i,p] <- a
      if(is.null(saic0))
      {
        saic0 <- a
        m0 <- tmp
      }
      else
      {
        if(a < saic0)
        {
          saic0 <- a
          m0 <- tmp
        }
      }
    }
  }
  if(length(dotsArgs)>0)
  {
    m1 <- carx(formula(m0), data=m0$data,p=m0$p,...)
    #outlier indices will be destroyed by the above command)
    m1$outlier.indices <- m0$outlier.indices
    m1$outlier.prefix <- m0$outlier.prefix
    m0 <- m1
  }
  list(fitted = m0, aicMat = aics)
}
