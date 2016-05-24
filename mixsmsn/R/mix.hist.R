##################################################################
##########     Plot histogram                #####################

mix.hist <- function(y, model, breaks=40, ...){
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Class of family",class(model),"not recognized.",sep=" "))
  y <- as.matrix(y) 
  if (dim(y)[2] != 1) stop("The mix.hist function is only appropriate for the univariate analysis.\n")
   
     if ((class(model) == "Skew.t") || class(model) == "t"){
      #### grafico ajuste
      hist(y, breaks = breaks,probability=T,col="grey",main=paste("Histogram of" , class(model),"fit"),...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col="red")
     }
     if (class(model) == "Skew.cn"){
      hist(y, breaks = breaks,probability=T,col="grey",main=paste("Histogram of" , class(model),"fit"),...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col="red")
     }
     if (class(model) == "Skew.slash"){
      hist(y, breaks = breaks,probability=T,col="grey",main=paste("Histogram of" , class(model),"fit"),...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col="red")
     }
     if ((class(model) == "Skew.normal") || (class(model) == "Normal")){
      hist(y, breaks = breaks,probability=T,col="grey",main=paste("Histogram of" , class(model),"fit"),...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape),col="red")
     }
}
