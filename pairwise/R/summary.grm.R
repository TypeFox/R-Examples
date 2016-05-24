#' @export summary.grm
#' @title S3 Summary for graphical Model Check
#' @description S3 summary method for object of class\code{"grm"}
#' @param object object of class\code{"grm"}
#' @param ci numeric with default \code{ci=2} to return cinfidence intervalls for point estimator.
#' @param ... other parameters passed trough
########################### hier die summary method #class grm #######################
summary.grm<-function(object,ci=2,...){
  #fu1<-function(x,ci){ cbind(x$parameter[,"sigma"]-(x$SE[,"sigma"]*ci) ,  x$parameter[,"sigma"]+(x$SE[,"sigma"]*ci))   }
  fu1<-function(x,ci){ cbind(x$sigma-(x$SEsigma*ci) ,  x$sigma+(x$SEsigma*ci))   }
  erg<-lapply(object,fu1,ci) 
  dimnames(erg[[1]])[[2]]<-paste(names(erg)[1] ,c("ci_l", "ci_u"))
  dimnames(erg[[2]])[[2]]<-paste(names(erg)[2] ,c("ci_l", "ci_u"))
  erg1<-cbind(erg[[1]],erg[[2]])
  u1<-erg1[,1]
  o1<-erg1[,2]
  u2<-erg1[,3]
  o2<-erg1[,4]
  
  OK <- ( (u1<u2 & o1>u2) | (u1<o2 & o1>o2) ) | ( (u1<u2 & o1>o2) | (u1>u2 & o1<o2) )
  
  #OK <- ((erg1[,3] > erg1[,1]) | (erg1[,4] > erg1[,1]) ) | ((erg1[,3] > erg1[,2]) | (erg1[,4] > erg1[,2]) )
  
  erg2<-cbind(erg1,OK)
  
  return(erg2)
}