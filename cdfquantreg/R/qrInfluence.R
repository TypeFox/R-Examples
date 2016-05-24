#' @title Influence Diagnosis For Fitted Cdfqr Object
#' @description Influence Diagnosis (dfbetas) For Fitted Cdfqr Object
#' @aliases influence.cdfqr
#' @param model A cdfqr model object 
#' @param type A string that indicates whether the results for all paramters are to be returned, or only the location/dispersion submodel's parameters returned.
#' @param what for influence statistics based on coefficient values, indicate the predictor variables that needs to be tested. 
#' @param plot if plot is needed.
#' @param method Currently only 'dfbeta' method is available.
#' @param ... currently ignored.s
#' @return A matrix, each row of which contains the estimated influence on parameters when that row's observation is removed from the sample.
#' 
#' @examples
#' data(cdfqrExampleData)
#'fit <- cdfquantreg(crc99 ~ vert | confl, 't2', 't2', data = JurorData)
#'influcne <- influence(fit)
#'plot(influcne[,2])
#'
#'\dontrun{
#' # Same as influence(fit)
#'dfbetval <- dfbetas(fit)
#'}
#'
#' @seealso \code{\link[stats]{lm.influence}}, \code{\link[stats]{influence.measures}}
#' @import stats
#' @method influence cdfqr
#' @import graphics
#' @export
influence.cdfqr <- function(model, method = "dfbeta",type = c("full","location", "dispersion"), 
                            what = "full",plot=FALSE, ...) {
  # - dfbeta residuals (for influence check)
 
  method <- match.arg(method)
  influcne <- switch(method, dfbeta = {
    infl.stats <- dfbeta(model, type, what)
  })
  
  return(influcne)
  if (plot){
    if (!is.null(ncol(influence))) {
      par(mfrow=c(3,2))
      for (i in 1:ncol(influence))
      {
        plot(influence[,i],xlab = colnames(influence)[i])
        text(influence[,i], labels=1:nrow(influence),cex=1, pos=2, col="red")
      }
    }
  }
}

#' @method dfbeta cdfqr
#' @export
#' @rdname influence.cdfqr
dfbeta.cdfqr <- function(model, type = c("full","location", "dispersion"), 
                          what = "full", ...) {
  dfbetas(model, type = type, what = what, ...)
  }

#' @method dfbetas cdfqr
#' @export
#' @rdname influence.cdfqr
dfbetas.cdfqr <- function(model, type = c("full","location", "dispersion"), what = "full", ...) {
  # - dfbeta residuals (for influence check)
  type <- match.arg(type)
  
  # The original model estimation call
  call <- model$call
  fd1 <- model$family$fd
  sd1 <- model$family$sd
  
  # Get the data from the call
  dat <- eval(call$data)
  n <- nrow(dat)
  
  # the number of location parameters
  k_lm <- nrow(model$coefficients$location)
  
  #Get the original coefficients estimates (mean and se)
  coefficients <- do.call(rbind, model$coefficients)
  coef0 <- coefficients[, 1]
  bse <- coefficients[, 2]
  
  # Casewise deletion
  betas <- NULL
  for (i in 1:nrow(dat)) {
    dat1 <- dat[-i, ]

    #modify the call for cdfqr function with the new dataset
    mod <- update(model, .~., fd=fd1, sd=sd1, data = dat1)
    
    dfbeta <- (coef(mod) - coef0)/bse
    betas <- rbind(betas, dfbeta)
  }
  
  colnames(betas) <- names(coef(model))
  betas_location <-betas[, 1:k_lm] 
  betas_dispersion <-betas[, (k_lm+1):ncol(betas)] 
  
  # If only a specific subset of parameters is needed, extract such subsets
  if (what != "full") {
    locind <- pmatch(colnames(betas_location), what, duplicates.ok = TRUE) 
    
    if(length(na.omit(locind))!=0) {
      betas_location<- betas_location[, what]
    }else{betas_location <- NULL}
  
    precind <- pmatch(colnames(betas_dispersion), what, duplicates.ok = TRUE) 
    
    if(length(na.omit(precind))!=0) {
      betas_dispersion<- betas_dispersion[, what]
    }else{betas_dispersion <- NULL}
    
  }
  
  #Rename the two submodels' paraters to distinguish the two submodels
  colnames(betas_location) <- paste("lm.",colnames(betas_location),sep="")
  colnames(betas_dispersion) <- paste("pm.",colnames(betas_dispersion),sep="")
  
  betas<-cbind(betas_location, betas_dispersion)
  
  betas <- switch(type, full = {
    betas
  }, location = {
    betas_location
  },  dispersion = {
    betas_dispersion
  })
  
  return(betas)
  
} 
