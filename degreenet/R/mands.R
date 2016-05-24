#  File degreenet/R/mands.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Calculate the modified Anderson-Darling Statistic
#
mands <- function(x,type="dp",v,cutoff=1,cutabove=1000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 n <- length(x)
 maxk <- max(x)
 mx <- cutoff:maxk
 cprob <- 1
#
# first calculates the ECDF
#
 ecdf <- cumsum(tabulate(x,nbins=maxk))/length(x)
 if(cutoff>1){
  ecdf <- ecdf[-c(1:(cutoff-1))]
 }
#
# now the theoretical distribution
#
 if(type=="dp"){
  cprob <- 1
  if(cutoff>1){
   cprob <- 1-sum(ddp(v,x=1:(cutoff-1)))
  }
  pdf <- ddp(v,x=mx) / cprob
  cdf <- cumsum(pdf)
 }
 if(type=="geom"){
  pdf <- dgeom(x=mx-cutoff, prob=1/v[1])
  cdf <- cumsum(pdf)
 }
 if(type=="nb"){
  pdf <- dnbinom(x=mx-cutoff, size=v[1]*v[2], prob=v[2])
  cdf <- cumsum(pdf)
 }
 if(type=="yule"){
  if(cutoff>1){
   cprob <- 1-sum(dyule(v,x=1:(cutoff-1)))
  }else{
   cprob <- 1
  }
  if(cutabove<1000){
   xc <- 1:cutabove
   aaa <- dyule(v,x=xc)
   cprob <- cprob - 1 + sum(aaa)
  }
  pdf <- dyule(v,x=mx) / cprob
  cdf <- cumsum(pdf)
 }
 if(type=="gy"){
  pdf <- dgyule(v,x=mx,cutoff=cutoff)
  cdf <- cumsum(pdf)
 }
 if(type=="war"){
  if(cutoff>1){
   cprob <- 1 - sum(dwar(v,x=1:(cutoff-1)))
  }else{
   cprob <- 1
  }
  pdf <- dwar(v,x=mx)/cprob
  cdf <- cumsum(pdf)
 }
 if(type=="pln"){
  if(cutoff>1){
   cprob <- 1 - sum(dpln(v,x=1:(cutoff-1)))
  }else{
   cprob <- 1
  }
  pdf <- dpln(v,x=mx)/cprob
  cdf <- cumsum(pdf)
 }
#if(type=="gpln"){
# pdf <- dgpln(v,x=mx,cutoff=cutoff)
# cdf <- cumsum(pdf)
#}
 if(type=="gwar"){
  pdf <- dgwar(v,x=mx,cutoff=cutoff)
  cdf <- cumsum(pdf)
 }
 if(type=="nbwar"){
  pdf <- dnbwar(v,x=mx,cutoff=cutoff)
  cdf <- cumsum(pdf)
 }
#
# funky distributions
#
 if(type=="pe"){
  xr <- 1:10000
  xr <- xr[xr >= cutoff]
  c0 <- sum(xr^(-v[1])*exp(-xr/v[2]))
  pdf <- exp(-mx/v[2])/(c0*mx^v[1])
  cdf <- cumsum(pdf)
 }
 if(type=="mix"){
  if(cutoff==1){
   aaa <- zeta(v[1])
  }else{
   aaa <- zeta(v[1])-sum((1:trunc(cutoff-0.999))^(-v[1]))
  }
  ppow <- 1/(aaa*mx^v[1])
  ppoi <- dpois(x=mx, lambda=v[2])/(ppois(cutoff-1,
lambda=v[2],lower.tail=FALSE))
  pdf <- v[3]*ppow + (1-v[3])*ppoi
  cdf <- cumsum(pdf)
 }
 if(type=="mixnb"){
  if(cutoff==1){
   aaa <- zeta(v[1])
  }else{
   aaa <- zeta(v[1])-sum((1:trunc(cutoff-0.999))^(-v[1]))
  }
  ppow <- 1/(aaa*mx^v[1])
  pdf <- dnbinom(x=mx, size=v[3], prob=v[2]/(1+v[2]))
  pdf <- pdf / (pnbinom(q=cutoff-1, size=v[3],
prob=v[2]/(1+v[2]),lower.tail=FALSE))
  pdf <- v[4]*ppow + (1-v[4])*pdf
  cdf <- cumsum(pdf)
 }
 if(type=="mixnbpoi"){
  ppoi <- dpois(x=mx, lambda=v[3])/(ppois(cutoff-1,
lambda=v[3],lower.tail=FALSE))
  pdf <- dnbinom(x=mx, size=v[2], prob=v[1]/(1+v[1]))
  pdf <- pdf / (pnbinom(q=cutoff-1, size=v[2],
prob=v[1]/(1+v[1]),lower.tail=FALSE))
  pdf <- v[4]*pdf + (1-v[4])*ppoi
  cdf <- cumsum(pdf)
 }
#
 if(hellinger){
#
# NOTE the next line
# Make the distribution the marginal and not conditional
#
  pdf <- cprob * pdf
#
  tr <- 0:max(x)
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff & tr <= cutabove]
# pc <- pc / sum(pc)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  ads <- -(out - 0.5 * (1-sum(pdf[pc>0]))) # penalized Hellinger Distance
 }else{
  ads <- n * sum( ( (ecdf-cdf)^2*pdf )/( cdf*(1-cdf) ),na.rm=TRUE)
 }
 eout <- cbind(mx,ecdf,cdf)
 dimnames(eout) <- list(NULL,c("k","ecdf","cdf"))
 list(cdf=eout, ads=ads)
}
#
# Calculate the modified Anderson-Darling Statistic
#
plotcdf <- function(v, maxk=15,type="dp",cutoff=1, ptype="p",add=FALSE){
 x <- x[x >= cutoff]
 n <- length(x)
 mx <- 1:maxk
 if(type=="dp"){
  if(cutoff==1){
   aaa <- zeta(v)
  }else{
   aaa <- zeta(v)-sum((1:trunc(cutoff-0.999))^(-v))
  }
  pdf <- 1/(aaa*mx^v)
  cdf <- cumsum(pdf)
 }
 if(type=="nb"){
  pdf <- dnbinom(x=mx, size=v[1]*v[2], prob=v[2])
  pdf <- pdf / (pnbinom(q=cutoff-1, size=v[2], prob=v[2],lower.tail=FALSE))
# pdf <- pdf/sum(pdf)
  cdf <- cumsum(pdf)
 }
 if(type=="mix"){
  if(cutoff==1){
   aaa <- zeta(v[1])
  }else{
   aaa <- zeta(v[1])-sum((1:trunc(cutoff-0.999))^(-v[1]))
  }
  ppow <- 1/(aaa*mx^v[1])
  ppoi <- dpois(x=mx, lambda=v[2])/(ppois(cutoff-1,
lambda=v[2],lower.tail=FALSE))
  pdf <- v[3]*ppow + (1-v[3])*ppoi
  cdf <- cumsum(pdf)
 }
 if(type=="mixnb"){
  if(cutoff==1){
   aaa <- zeta(v[1])
  }else{
   aaa <- zeta(v[1])-sum((1:trunc(cutoff-0.999))^(-v[1]))
  }
  ppow <- 1/(aaa*mx^v[1])
  pdf <- dnbinom(x=mx, size=v[3], prob=v[2]/(1+v[2]))
  pdf <- pdf / (pnbinom(q=cutoff-1, size=v[3],
prob=v[2]/(1+v[2]),lower.tail=FALSE))
  pdf <- v[4]*ppow + (1-v[4])*pdf
  cdf <- cumsum(pdf)
 }
 if(type=="mixnbpoi"){
  ppoi <- dpois(x=mx, lambda=v[3])/(ppois(cutoff-1,
lambda=v[3],lower.tail=FALSE))
  pdf <- dnbinom(x=mx, size=v[2], prob=v[1]/(1+v[1]))
  pdf <- pdf / (pnbinom(q=cutoff-1, size=v[2],
prob=v[1]/(1+v[1]),lower.tail=FALSE))
  pdf <- v[4]*pdf + (1-v[4])*ppoi
  cdf <- cumsum(pdf)
 }
 if(type=="mixturenbpoi"){
  plog <- ""
  ptype <- "l"
  if(v[4] > 0.999){mmx <- qgamma(0.999,shape=v[2], scale=v[1])}
  if(v[4] < 0.1){mmx <- v[3]+1}
  if(0.1 < v[4] & v[4] < 0.999){mmx <- max(v[3]+1, qgamma(0.999,shape=v[2], scale=v[1]))}
  mx <- c(v[3],seq(0.005,mmx,length=100))
  pdf <- dgamma(x=mx, shape=v[2], scale=v[1])
 print(c(v,pgamma(1,shape=v[2], scale=v[1])))
  pdf <- v[4]*pdf
  pdf[1] <- pdf[1] - (1-v[4])/diff(mx)[1]
  pdf <- pdf[order(mx)]
  mx <- mx[order(mx)]
  cdf <- cumsum(pdf)*diff(-mx)
 }
 if(type=="yule"){
  pdf <- exp(log(v-1) + lgamma(mx) + lgamma(v) - lgamma(mx+v) )
  cdf <- cumsum(pdf)
 }
 if(type=="pe"){
  xr <- 1:10000
  xr <- xr[xr >= cutoff]
  aaa <- sum(xr^(-v[1])*exp(-xr/v[2]))
  pdf <- exp(-mx/v[2])/(aaa*mx^v[1])
  cdf <- cumsum(pdf)
 }
#
 if(add){
  if(ptype=="l"){
   lines(x=mx,y=pdf,xlab="degree",ylab="density",log=plog)
  }else{
   points(x=mx,y=pdf,xlab="degree",ylab="density",log=plog)
  }
 }else{
  plot(x=mx,y=pdf,xlab="degree",ylab="density",type=ptype,
       log=plog)
 }
 cbind(mx,pdf,cdf)
}
#
# MANDS bootstrap functions
#
bsyule <- function(x, cutoff=1, m=200, np=1, alpha=0.95, v=NULL,
                   hellinger=FALSE, cutabove=1000){
 if(missing(v)){v <- ayulemle(x=x,cutoff=cutoff,cutabove=cutabove,
                      hellinger=hellinger)$theta}
 obsmands <- mands(x,type="yule",v=v,cutoff=cutoff,cutabove=cutabove,
                   hellinger=hellinger)
 bsmles <- matrix(0,nrow=m,ncol=np+1)
 for(i in 1:m){
#
# Sample from a Yule CDF at the MLE
#
  xsamp <- sample(x=obsmands$cdf[,"k"], 
                  size=length(x),
                  replace=TRUE,
		  prob=diff(c(0,obsmands$cdf[,"cdf"]))
		 )
  bsmles[i,-np-1] <- ayulemle(x=xsamp,cutoff=cutoff,cutabove=cutabove,
                      hellinger=hellinger)$theta
  bsmles[i, np+1] <- mands(xsamp,type="yule",v=bsmles[i,-np-1],
                      cutoff=cutoff,cutabove=cutabove,
                      hellinger=hellinger)$ads
 }
 if(hellinger){
  dimnames(bsmles) <- list(1:m,c("PDF MLE","Pen. Hell. Dist."))
 }else{
  dimnames(bsmles) <- list(1:m,c("PDF MLE","MANDS"))
 }
 #
 list(dist=obsmands$cdf,
  obsmle=v,
  bsmles=bsmles,
  quantiles=quantile(bsmles[,np+1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
  pvalue=sum(obsmands$ads < bsmles[,np+1],na.rm=TRUE)/(sum(bsmles[,np+1]>0,na.rm=TRUE)+1),
  obsmands=obsmands$ads,
  meanmles=apply(bsmles,2,mean)
 )
}
bswar <- function(x, cutoff=1, m=200, np=2, alpha=0.95, v=NULL,
                   hellinger=FALSE){
 if(missing(v)){v <- awarmle(x=x,cutoff=cutoff,hellinger=hellinger)$theta}
 obsmands <- mands(x,type="war",v=v,cutoff=cutoff,hellinger=hellinger)
 bsmles <- matrix(0,nrow=m,ncol=np+1)
 for(i in 1:m){
#
# Sample from a Yule CDF at the MLE
#
  xsamp <- sample(x=obsmands$cdf[,"k"], 
                  size=length(x),
                  replace=TRUE,
		  prob=diff(c(0,obsmands$cdf[,"cdf"]))
		 )
  wmle <- awarmle(x=xsamp,cutoff=cutoff,hellinger=hellinger)$theta
  bsmles[i,-np-1] <- wmle
  bsmles[i, np+1] <- mands(xsamp,type="war",v=wmle,
                      cutoff=cutoff,hellinger=hellinger)$ads
# if(done){print(paste("done",i))}
 }
 if(hellinger){
  dimnames(bsmles) <- list(1:m,c("PDF MLE","prob. new","Pen. Hell. Dist."))
 }else{
  dimnames(bsmles) <- list(1:m,c("PDF MLE","prob. new","MANDS"))
 }
 #
 list(dist=obsmands$cdf,
  obsmle=v,
  bsmles=bsmles,
  quantiles=quantile(bsmles[,np+1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
  pvalue=sum(obsmands$ads < bsmles[,np+1],na.rm=TRUE)/(sum(bsmles[,np+1]>0,na.rm=TRUE)+1),
  obsmands=obsmands$ads,
  meanmles=apply(bsmles,2,mean)
 )
}
bsdp <- function(x, cutoff=1, m=200, np=1, alpha=0.95){
 mle <- adpmle(x=x,cutoff=cutoff)$theta
 obsmands <- mands(x,type="dp",v=mle,cutoff=cutoff)
 bsmles <- matrix(0,nrow=m,ncol=np+1)
 for(i in 1:m){
#
# Sample from a Power CDF at the MLE
#
  xsamp <- sample(x=obsmands$cdf[,"k"], 
                  size=length(x),
                  replace=TRUE,
		  prob=diff(c(0,obsmands$cdf[,"cdf"]))
		 )
  bsmles[i,-np-1] <- adpmle(x=xsamp,cutoff=cutoff)$theta
  bsmles[i, np+1] <- mands(x=xsamp,type="dp",v=bsmles[i,-np-1],cutoff=cutoff)$ads
 }
 dimnames(bsmles) <- list(1:m,c("PDF MLE","MANDS"))
 #
 list(dist=obsmands$cdf,
  obsmle=mle,
  bsmles=bsmles,
  quantiles=quantile(bsmles[,np+1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
  pvalue=sum(obsmands$ads < bsmles[,np+1],na.rm=TRUE)/(sum(bsmles[,np+1]>0,na.rm=TRUE)+1),
  obsmands=obsmands$ads,
  meanmles=apply(bsmles,2,mean)
 )
}
#
bsnb <- function(x, cutoff=1, m=200, np=2, alpha=0.95,
                   hellinger=FALSE){
 mle <- anbmle(x=x,cutoff=cutoff,hellinger=hellinger)$theta
 obsmands <- mands(x,type="nb",v=mle,cutoff=cutoff,hellinger=hellinger)
 bsmles <- matrix(0,nrow=m,ncol=np+1)
 for(i in 1:m){
#
# Sample from a Yule CDF at the MLE
#
  xsamp <- sample(x=obsmands$cdf[,"k"], 
                  size=length(x),
                  replace=TRUE,
		  prob=diff(c(0,obsmands$cdf[,"cdf"]))
		 )
  bsmles[i,-np-1] <- anbmle(x=xsamp,cutoff=cutoff,hellinger=hellinger)$theta
  bsmles[i, np+1] <- mands(xsamp,type="nb",v=bsmles[i,-np-1],cutoff=cutoff,hellinger=hellinger)$ads
 }
 if(hellinger){
  dimnames(bsmles) <- list(1:m,c("e.stop MLE","pr.1stop","Pen. Hell. Dist."))
 }else{
  dimnames(bsmles) <- list(1:m,c("expected count","Prob. of a stop","MANDS"))
 }
 #
 list(dist=obsmands$cdf,
  obsmle=mle,
  bsmles=bsmles,
  quantiles=quantile(bsmles[,np+1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
  pvalue=sum(obsmands$ads < bsmles[,np+1],na.rm=TRUE)/(sum(bsmles[,np+1]>0,na.rm=TRUE)+1),
  obsmands=obsmands$ads,
  meanmles=apply(bsmles,2,mean)
 )
}
