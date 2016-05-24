findLogMargLik <-
function (model, tData, dimens, alpha) {
  
  model <- na.omit(model)
  logMargLik <- 0  

  .C("findLogMargLik", as.integer(tData[model,]), as.integer(dim(tData)[2]), as.integer(length(model)), as.integer(dimens[model]), as.double(alpha), logMargLik = as.double(logMargLik), package = "c_code")$logMargLik
  
}
