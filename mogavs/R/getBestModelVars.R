getBestModelVars <-
function(mogavs,nvars,data,method=c("AIC","BIC","mse",NULL)){
  return(as.vector(sapply(nvars, function(x)colnames(data)[which(getBestModel(mogavs,x,method=method)==1)+1])))
}
