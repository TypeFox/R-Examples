##################################################################
##########     Plot density                  #####################
mix.dens <- function(y, model, log=FALSE, ylab=NULL, xlab = NULL, main = NULL, ...){
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Class of family",class(model),"not recognized.",sep=" "))
  y <- as.matrix(y) 
  if (dim(y)[2] != 1) stop("The mix.dens function is only appropriate for the univariate analysis.\n")

     if(length(ylab) == 0) ylab = class(model)
     if(length(xlab) == 0) xlab = "x"

     if(length(main) == 0){
       main <- "density plot"
       if(log) main <- "log-density plot"
     }
     
     lim <- -9

     xx=seq(min(y),max(y),(max(y)-min(y))/1000)   
     if ((class(model) == "Skew.t") || class(model) == "t"){
      if(!log) plot(xx,d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),type="l",ylab=ylab, xlab = xlab, main = main, ...)
      else{
         aux <- log(d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         plot(xx,aux,type="l",ylab=ylab, xlab = xlab, main = main, ...)
	}
     }
     if (class(model) == "Skew.cn"){
      if(!log) plot(xx,d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),type="l",ylab=ylab, xlab = xlab, main = main, ...)
      else{
	 aux <- log(d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         plot(xx,aux,type="l",ylab=ylab, xlab = xlab, main = main, ...)
	}
     }
     if (class(model) == "Skew.slash"){
      if(!log) plot(xx,d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),type="l",ylab=ylab, xlab = xlab, main = main, ...)
      else{
	 aux <- log(d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         plot(xx,aux,type="l",ylab=ylab, xlab = xlab, main = main, ...)
	}
     }
     if ((class(model) == "Skew.normal") || (class(model) == "Normal")){
      if(!log) plot(xx,d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape),type="l",ylab=ylab, xlab = xlab, main = main, ...)
      else{
	 aux <- log(d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape))
	 aux[which(aux < lim)] <- NA
         plot(xx,aux,type="l",ylab=ylab, xlab = xlab, main = main, ...)
	}
     }
}
##################################################################
##########     Lines	                     #####################

mix.lines <- function(y, model, log=FALSE, ...){
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Class of family",class(model),"not recognized.",sep=" "))
  y <- as.matrix(y) 
  if (dim(y)[2] != 1) stop("The mix.dens function is only appropriate for the univariate analysis.\n")

     lim <- -9
     xx=seq(min(y),max(y),(max(y)-min(y))/1000)   

     if ((class(model) == "Skew.t") || class(model) == "t"){
      if(!log) lines(xx,d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu), ...)
      else{
         aux <- log(d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         lines(xx,aux,...)
	}
     }
     if (class(model) == "Skew.cn"){
      if(!log) lines(xx,d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),...)
      else{
	 aux <- log(d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         lines(xx,aux,...)
	}
     }
     if (class(model) == "Skew.slash"){
      if(!log) lines(xx,d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),...)
      else{
	 aux <- log(d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu))
	 aux[which(aux < lim)] <- NA
         lines(xx,aux,...)
	}
     }
     if ((class(model) == "Skew.normal") || (class(model) == "Normal")){
      if(!log) lines(xx,d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape),...)
      else{
	 aux <- log(d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape))
	 aux[which(aux < lim)] <- NA
         lines(xx,aux,...)
	}
     }
}
