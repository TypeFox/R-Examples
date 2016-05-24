if(getRversion() >= "2.15.1") utils::globalVariables("loess.sd")


g.dl6<-function(x,l, p1, p2){
	dlvalue <- rep(0.0, length(l))
	res <- .C("doDL", x, l, p1, p2, length(l), dl_values=dlvalue)
	return(res$dl_values)
}

loess.lookup<-function(look){
   #-- spline function
#     a<- 2.1828214 
#     b<- -0.0335983
#     c<- 0.0463050 
#     d<- -0.0495180 
#   
#     x1<- 51.2058
#     x2<-65.801 
#   
#     val<-a+b*look+c*(ifelse(look>x1,look-x1,0))+d*(ifelse(look>x2,look-x2,0))
#     
#     if(look>77.2) val<-a+b*77.2+c*(ifelse(77.2>x1,77.2-x1,0))+d*(ifelse(77.2>x2,77.2-x2,0))
      
   # call data(loess_sd) before using this function
   idx <- cut(look, loess.sd$x, labels=FALSE, include.lowest = TRUE)
   return(loess.sd$y[idx])
}


dnorm.trunc<-function(x,mean,sd,low,high){
  out<-dnorm(x,mean=mean,sd=sd)/(pnorm(high,mean=mean,sd=sd)-pnorm(low,mean=mean,sd=sd))
  out[x<low]<-0
  out[x>high]<-0
  return(out)
}

rnorm.trunc<-function(mean,sd,low,high){
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low || temp>high) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else high
  	warning(paste('Maximum iterations reached in rnorm.trunc(', 
  				mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
  }
  return(temp)
}

rgamma.ltrunc<-function(shape,rate,low){
  temp<--999
  maxit <- 10
  i <- 1
  while(temp<low && i <= maxit) {
     temp<-rgamma(1,shape=shape,rate=rate)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- low
  	warning(paste('Maximum iterations reached in rgamma.ltrunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}


