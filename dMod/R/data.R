#' Compare data and model prediction by computing residuals
#' 
#' @param data data.frame with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param out output of ode(), optionally augmented with attributes 
#' "deriv" (output of ode() for the sensitivity equations) and
#' "parameters" (character vector of parameter names, a subsest of those 
#' contained in the sensitivity equations). If "deriv" is given, also "parameters"
#' needs to be given.
#' @return data.frame with the original data augmented by columns "prediction" (
#' numeric, the model prediction), "residual" (numeric, difference between
#' prediction and data value), "weighted.residual" (numeric, residual devided
#' by sigma). If "deriv" was given, the returned data.frame has an 
#' attribute "deriv" (data.frame with the derivatives of the residuals with 
#' respect to the parameters).
#' @export
#' @import cOde
res <- function (data, out) {
  
  # Unique times, names and parameter names
  times <- sort(unique(data$time))
  names <- as.character(unique(data$name))
  pars <- attr(out, "parameters")
  
  # Match data times/names in unique times/names
  data.time <- match(data$time, times)
  data.name <- match(data$name, names)
  
  # Match unique times/names in out times/names
  time.out <- match(times, out[,1])
  name.out <- match(names, colnames(out))
  
  # Match data times/names in out times/names
  timeIndex <- time.out[data.time]
  nameIndex <- name.out[data.name]
  prediction <- sapply(1:nrow(data), function(i) out[timeIndex[i], nameIndex[i]]) 
  
  # Propagate derivatives if available
  deriv <- attr(out, "deriv")
  deriv.data <- NULL
  if (!is.null(deriv)) {
    sensnames <- as.vector(outer(names, pars, paste, sep="."))
    # Match names to the corresponding sensitivities in sensnames
    names.sensnames <- apply(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)), 1, identity)
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv))
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
    # Derivatives of the prediction
    deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) deriv[timeIndex[i], derivnameIndex[, i]]))
    colnames(deriv.prediction) <- pars
    
    deriv.data <- data.frame(time = data$time, name = data$name, deriv.prediction)
  }
  
  # Compute residuals
  residuals <- prediction - data$value 
  weighted.residuals <- (prediction - data$value)/data$sigma
  data <- cbind(data, prediction = prediction, residual = residuals, 
                weighted.residual = weighted.residuals)
  data <- data[c("time", "name", "value", "prediction", "sigma", 
                 "residual", "weighted.residual")]
  attr(data, "deriv") <- deriv.data
  return(data)
}


#' Compute the weighted residual sum of squares
#' 
#' @param nout data.frame (result of \link{res})
#' @return list with entries value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric).
#' @export
wrss <- function(nout) {
  
  obj <- sum(nout$weighted.residual^2)
  grad <- NULL
  hessian <- NULL
  
  if(!is.null(attr(nout, "deriv"))) {
    nout$sigma[is.na(nout$sigma)] <- 1 #replace by neutral element
  
    sens <- as.matrix(attr(nout, "deriv")[,-(1:2)])
    grad <- as.vector(2*matrix(nout$residual/nout$sigma^2, nrow=1)%*%sens)
    names(grad) <- colnames(sens)
    hessian <- 2*t(sens/nout$sigma)%*%(sens/nout$sigma)
    
    
  }
  
  out <- list(value=obj, gradient=grad, hessian=hessian)
  class(out) <- c("obj", "list")
  
  return(out)

}


#' Generate dummy list of class \code{obj} from named numeric
#' 
#' @param p Names numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.obj(p)
#' @export
as.obj <- function(p) {
  
  obj <- list(
    value = 0,
    gradient = structure(rep(0, length(p)), names = names(p)),
    hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
  class(obj) <- "obj"
  
  return(obj)
  
}

