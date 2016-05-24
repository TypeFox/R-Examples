
param <- function(x,...) UseMethod("param") 


param.dyntermstrc_nss <- function(x,...) {
  param <- list()
  for(i in seq(x[[1]]$n_group)) param[[i]] =  t(mapply(function(j) x[[j]]$opt_result[[i]]$par,seq_along(x)))
  names(param) <- x[[1]]$group                          
  class(param) <- "dyntermstrc_param"
  param
}


param.dyntermstrc_yields <- function(x,...){
    param <- list() 
    param[[1]] <- x$optparam
    class(param) <- "dyntermstrc_param"
    param
}

summary.dyntermstrc_param <- function(object,type="none",lags=1,selectlags="Fixed", ...) {
    x <- object
    sumry <- list()
    length(sumry) <- length(x) 
    for(i in seq_along(x)) {
    
    # Augmented Dickey Fuller Test for levels
    sumry[[i]]$adflevels <- apply(x[[i]],2,function(x) ur.df(x,type=type,lags=lags,selectlags=selectlags)) #alternatively use adf.testx  

    sumry[[i]]$adflevelsm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x[[i]]))

    for (j in 1:length(sumry[[i]]$adflevels)) {
      sumry[[i]]$adflevelsm[switch(type,"none"=1,"trend"=c(1:3)),j] <- sumry[[i]]$adflevels[[j]]@teststat # adf.test : $statisic
      sumry[[i]]$adflevelsm[switch(type,"none"=2,"trend"=4),j] <- sumry[[i]]$adflevels[[j]]@lags # adf.test: $parameter
      sumry[[i]]$adflevelsm[switch(type,"none"=3,"trend"=c(5:7)),j] <- sumry[[i]]$adflevels[[j]]@cval[,3] # adf.test: $p.value
    }
    rownames(sumry[[i]]$adflevelsm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
    colnames(sumry[[i]]$adflevelsm) <- colnames(x[[i]])
  
  
    # Augmented Dickey Fuller Test for first differences
    sumry[[i]]$adfdiff <- apply(x[[i]],2,function(x) ur.df(diff(x),type=type,lags=lags,selectlags=selectlags))
    sumry[[i]]$adfdiffm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x[[i]]))
    for (j in 1:length(sumry[[i]]$adflevels)) {
      sumry[[i]]$adfdiffm[switch(type,"none"=1,"trend"=c(1:3)),j] <- sumry[[i]]$adfdiff[[j]]@teststat # adf.test : $statisic
      sumry[[i]]$adfdiffm[switch(type,"none"=2,"trend"=4),j] <- sumry[[i]]$adfdiff[[j]]@lags # adf.test: $parameter
      sumry[[i]]$adfdiffm[switch(type,"none"=3,"trend"=c(5:7)),j] <- sumry[[i]]$adfdiff[[j]]@cval[,3] # adf.test: $p.value
    }
    rownames(sumry[[i]]$adfdiffm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
    colnames(sumry[[i]]$adfdiffm) <- colnames(x[[i]])

    sumry[[i]]$paramcor <- cor(x[[i]])
    sumry[[i]]$diffparamcor <- cor(apply(x[[i]],2,diff))

  }
  names(sumry) <- names(x)   
    
  class(sumry) <- "summary.dyntermstrc_param"
  sumry
}


print.summary.dyntermstrc_param <- function(x, ...) {
  for(i in seq_along(x)) {
  cat("---------------------------------------------------\n")
  cat(paste("ADF for ",names(x)[[i]],": ",sep=""))
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  # Augmented Dickey Fuller Test for levels
  print.default(t(x[[i]]$adflevelsm))
  cat("\n")
  cat("---------------------------------------------------\n")
  cat(paste("ADF of differences for ",names(x)[[i]],": ",sep=""))
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  # Augmented Dickey Fuller Test for first differences
  print.default(t(x[[i]]$adfdiffm))
  cat("\n")
  # correlation matrix of parameters
  cat("---------------------------------------------------\n")
  cat(paste("Correlation of parameters for ",names(x)[[i]],": ",sep=""))
  cat("\n")
  cat("---------------------------------------------------\n")
  print.default(x[[i]]$paramcor)
  cat("\n")
  cat("---------------------------------------------------\n")
  cat(paste("Correlation of differences for ",names(x)[[i]],": ",sep=""))
  cat("\n")
  cat("---------------------------------------------------\n")
  print.default(x[[i]]$diffparamcor)
  cat("\n")
}
}


plot.dyntermstrc_param <- function(x,type="param",...){
  old.par <- par(no.readonly = TRUE) 
  
  
  # 2D plot of parameters
  if(type=="param") {
    if(ncol(x[[1]])<=3) mfrow = c(1,3)
    if(ncol(x[[1]])==4) mfrow = c(2,2)
    if(ncol(x[[1]])>4 && ncol(x[[1]]) <= 6) mfrow = c(2,3)

    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)
        
   for(i in seq_along(x)){
  
    param <- x[[i]]
    
   
    plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
                col=1,lwd=2,... )
           grid()
           plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
           grid()
           plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
           grid()
    
    if(ncol(param)==4) {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
    }
    
    if(ncol(param)==6) {
           plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
           col=4,lwd=2,... )
           grid()
           plot(param[,5],type="l",xlab="Time",ylab=expression(hat(beta)[3]),
           col=5,lwd=2,... )
           grid()
           plot(param[,6],type="l",xlab="Time",ylab=expression(hat(tau)[2]),
           col=6,lwd=2,... )
           grid()
    }
   } 
  }

 # 2D plot of parameter differences
  if(type=="diffparam") {
     if(ncol(x[[1]])<=3) mfrow = c(1,3)
    if(ncol(x[[1]])==4) mfrow = c(2,2)
    if(ncol(x[[1]])>4 && ncol(x[[1]]) <= 6) mfrow = c(2,3)

    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)

    for(i in seq_along(x)){
      param <- x[[i]]
      diffparam <- apply(param,2,diff)

      for(i in seq(ncol(diffparam))) {
       
        plot(diffparam[,i],type="l",xlab="Time",
           ylab= paste("delta",colnames(diffparam)[i],sep=" "),col=i,lwd=2,... )
        grid()
      }
    }
  }

    # ACF/PCF
  if(type=="acf") {
    if(ncol(x[[1]])==3) mfrow = c(2,3)
    if(ncol(x[[1]])==4) mfrow = c(4,2)
    if(ncol(x[[1]])==6) mfrow = c(4,3)

    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)
    for(i in seq_along(x)){
      param <- x[[i]]
      if(ncol(param) > 3 ){
       for(i in 1:(ncol(param)/2)) acf(param[,i],main=colnames(param)[i])
       for(i in 1:(ncol(param)/2)) pacf(param[,i],main=colnames(param)[i])
    
       for(i in (ncol(param)/2+ 1):ncol(param)) acf(param[,i],main=colnames(param)[i])
       for(i in (ncol(param)/2+ 1):ncol(param)) pacf(param[,i],main=colnames(param)[i])
      } else {

       for(i in 1:ncol(param)) acf(param[,i],main=colnames(param)[i])
       for(i in 1:ncol(param)) pacf(param[,i],main=colnames(param)[i])

      }
    } 
  }

   
  on.exit(par(old.par))


}

fcontrib <- function(x, method="ns",lambda=0.0609*12, index=1, m=1:10, ylim=NULL ,... ) UseMethod("fcontrib")

fcontrib.dyntermstrc_param <- function(x, method="ns",lambda=0.0609*12, index=1, m=1:10, ylim=NULL ,... ){
  par(if(length(x) > 1) mfrow=length(x),... )
  for(i in seq_along(x)){
    param <- x[[i]]
    if(ncol(param)==3) param <- cbind(param,1/lambda)
    fc1 <- param[,1]
    fc2 <- t(mapply(function(i) param[i,2]*((1-exp(-m/param[i,4]))/(m/param[i,4])), seq(nrow(param))))
    fc3 <- t(mapply(function(i) param[i,3]*(((1-exp(-m/param[i,4]))/(m/param[i,4]))-exp(-m/param[i,4])), seq(nrow(param))))
    if(ncol(param)==6) fc4 <- t(mapply(function(i) param[i,5]*(((1 - exp(-m/param[i,6]))/(m/param[i,6])) - exp(-m/param[i,6])), seq(nrow(param))))
    
    if(is.null(ylim)) ylim <- c(min(fc1,fc2,fc3),max(fc1,fc2,fc3))          
    
    plot(m,rep(fc1[index],length(m)), col=1,type="l",lty=1, ylim=ylim, xlab="Time to maturity", ylab="Factor contribution",lwd=2,main=get_realnames(method))
    
    
    ## beta_1*( )
    lines(m, fc2[index,],type="l",col=2,lty=3,lwd=2)
    ## beta_2*()
    lines(m, fc3[index,],lty=4,col=4,lwd=2)
    ## beta_3*()
    if(ncol(param)==6) lines(m, fc4[index,],lty=5,col=5,lwd=2)
    
    abline(h=0,lty=1, lwd = 1, col = "grey")
     
    legend("topright",bg='white', box.col="white", box.lwd = 0,
           legend=c(expression(beta[0]),
             expression(beta[1]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1])))),
             expression(beta[2]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1]))-exp(-frac(m,tau[1])))),
             if(method=="sv") expression(beta[3]*(frac(1-exp(-frac(m,tau[2])),frac(m,tau[2]))-exp(-frac(m,tau[2])))),
             if(method=="asv")  expression(beta[3]*(frac(1-exp(-frac(m,tau[2])),frac(m,tau[2]))-exp(-frac(2*m,tau[2]))))                 
             ),
           lty=c(1,2,4,5), col=c(1,2,4,5),bty="o", cex = 0.9
           ) 
    box()  
  }
}
