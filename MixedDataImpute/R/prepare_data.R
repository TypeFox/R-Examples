#' Prepare a dataset for imputation
#'
#' Prepares a dataset for imputation by mapping factor levels to integers and scaling
#' Y. Primarily used by \code{hcmm_impute} internally
#'
#' @param X A data frame of categorical variables (as factors)
#' @param Y A matrix or data frame of continuous variables
#' @param init If \code{TRUE}, initialize missing values (from a marginal bootstrap of
#' observed values)
#'
#' @return An object of class \code{hcmm_data}
#' @export
#'
prepare_data = function(X, Y, init=TRUE) {
  
  if(nrow(X)!=nrow(Y)) stop("X and Y must have the same number of observations")
  if(!all(sapply(X, function(tt) is.factor(tt))))  stop("X must contain factors")
  if(!all(sapply(Y, function(tt) is.numeric(tt)))) stop("Y must have numeric columns")
  
  Xmap = lapply(X, function(x) mapLevels(x))
  cx = sapply(Xmap, length)
  Xint = data.frame(lapply(X, as.integer))
  Yscaled = scale(Y)
  y.m = attr(Yscaled, 'scaled:center')
  y.s = attr(Yscaled, 'scaled:scale')
  
  if(init) {
    for(i in 1:ncol(Yscaled)) {
      Yscaled[is.na(Yscaled[,i]),i] = sample(Yscaled[!is.na(Yscaled[,i]),i], 
                                         sum(is.na(Yscaled[,i])),
                                         replace=TRUE)
    }
    
    for(i in 1:ncol(Xint)) {
      Xint[is.na(Xint[,i]),i] = sample(Xint[!is.na(Xint[,i]),i], 
                                             sum(is.na(Xint[,i])),
                                             replace=TRUE)
    }
  }
  
  ret = list(
    X=X, Y=Y,cx=cx,Xmap=Xmap,
    Xint=Xint-1, Yscaled=Yscaled,
    y.mean = y.m, y.scale = y.s
  )
  
  class(ret) = "hcmm_data"
  
  ret
  
}

#' Map raw imputations back to original scale/factor labels.
#'
#' Map raw imputations back to original scale (for continuous data) or factor labels. 
#' Most users can ignore this function, which is primarily used by \code{hcmm_impute} 
#' internally.
#'
#' @param Ximp A raw imputed X matrix
#' @param Yimp A raw imputed Y matrix
#' @param hcmmdat An \code{hcmm_data} object, used to recode/rescale the raw imputations
#'
#' @return A list with elements \code{X} and \code{Y} (transformed imputations).
#' @export
#'
remap_imputations = function(Ximp, Yimp, hcmmdat) {
  Xi = Ximp + 1
  Yi = t(hcmmdat$y.scale*t(Yimp)+hcmmdat$y.mean)
  Xn = as.data.frame(lapply(1:ncol(Xi), function(x) { xx = as.integer(Xi[,x]); mapLevels(xx)<-hcmmdat$Xmap[[x]]; xx }))
  
  colnames(Xn) = colnames(hcmmdat$X)
  colnames(Yi) = colnames(hcmmdat$Y)
  
  return(list(X=Xn, Y=Yi))
}