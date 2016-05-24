#' @export summary.S2CFA
#' @S3method summary S2CFA
#' @title S3 Summary for S2CFA
#' @description S3 summary method for object of class\code{"S2CFA"}
#' @param object object of class\code{"S2CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character with default \code{type="ex.fisher.test"}, to return wether the observed pattern are 'discriminating Types' or not significant at all based on the respective p-value. Another option for \code{type} is \code{type="pChi"}.     
#' @param ... other parameters passed trough

########################### hier die summary method #class S2CFA #######################
summary.S2CFA<-function(object, digits=3, type="ex.fisher.test",...){
  local.test <- object$local.test
  global.test <- object$global.test
  
  # object$bonferroni.alpha
  disc.Type <- ifelse(local.test[,which(names(local.test)==type)] < object$bonferroni.alpha,yes="+", no=".")
  
  templocal <- data.frame(pat.=local.test[,1], disc.Type, round(local.test[,2:3],digits=digits), round(local.test[,5:10],digits=digits) , check.names = FALSE)
  #print(templocal)
  # cat("results of global tests:","currently not implemented !","\n")
  cat("results of local tests: ","\n")
  cat("discriminating Type (+) / not discriminating Type (.) based on:", type, "; Bonferoni adj. alpha:", object$bonferroni.alpha,"\n")
  return(templocal)
}