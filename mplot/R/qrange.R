#' Identify an appropriate rance of values over which to bootstrap
#'
#' This function takes an upper and lower dimension size (obtained by
#' forwards and backwards model selection and then adding and subtracting
#' 2 from each of the extremes to encompas a broader range of models).
#' For both the small and large model size, the "best" model is identified
#' using the \code{leaps} package and the corresponding lack of fit measure
#' is calculated.
#'
#' @param k.range list with dimension elements k.max and k.min
#' @param yname name of dependent variable
#' @param fixed the full model formula
#' @param data full data table
#' @param method method used in Qm
#' @param force.in which variables to force into the model
#' @param model.type currently only lm or glm
#' @param family for glms.
#' @noRd
qrange = function(k.range,yname,fixed,
                  data,method,force.in,
                  model.type,family){
  kf = k.range$k.max
  if(model.type=="lm"){
    cand.models = summary(leaps::regsubsets(x = fixed,
                                     data = data,
                                     nbest = 1,
                                     nvmax = kf,
                                     force.in=force.in))$which+0
  } else if(model.type=="glm"){
    cand.models = bestglm::bestglm(Xy=data,
                          family=family,
                          IC = "BIC",
                          TopModels = 1,
                          nvmax = )$Subsets[,1:kf]+0
  }
  small.row = which(rowSums(cand.models)==k.range$k.min)
  small.names = colnames(cand.models)[which(cand.models[small.row,]==1)]
  if(length(small.names)>1){
    small.ff = paste(yname," ~ ",paste(small.names[-1],collapse="+"),sep="")
  } else small.ff = paste(yname,"~1")
  small.ff = stats::as.formula(small.ff)

  big.row = which(rowSums(cand.models)==kf)
  big.names = colnames(cand.models)[which(cand.models[big.row,]==1)]
  if(length(big.names)>1){
    big.ff = paste(yname," ~ ",paste(big.names[-1],collapse="+"),sep="")
  } else big.ff = paste(yname,"~1")
  big.ff = stats::as.formula(big.ff)
  if(model.type=="lm"){
    small.em = stats::lm(small.ff, data=data)
    big.em = stats::lm(big.ff, data=data)
  } else if(model.type=="glm"){
    small.em = stats::glm(small.ff, data=data, family=family)
    big.em = stats::glm(big.ff, data=data, family=family)
  }
  Q.min = Qm(big.em,method=method)
  Q.max = Qm(small.em,method=method)
  return(list(Q.min=Q.min,Q.max=Q.max))
}
