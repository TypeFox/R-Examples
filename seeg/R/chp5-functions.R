# empirical and theoretical normal cumulative plots

cdf.plot.gof <- function(x,dist="normal", mu=0,sd=1,rate=1,min=0,max=1,shape=1,scale=1){
  nx <- length(x)		# number of observations
  vx <- sort(x)		# data sorted
  Fx <- seq(nx)/nx	# ranks as fractions (quantiles)
  panel2(size=7)
  plot(vx, Fx, xlab="x", ylab="F(x)", ylim=c(0,1)) 	# empirical cumulative plot
  vnx <- seq(vx[1],vx[nx],length.out=100) # sequence with increased res
  # cdf values
  if(dist=="normal") {Fxh <- pnorm(vnx,mu,sd);Fxe <- pnorm(vx,mu,sd)}  		
  if(dist=="exp") {Fxh <- pexp(vnx,rate);Fxe <- pexp(vx,rate)}
  if(dist=="unif") {Fxh <- punif(vnx,min,max);Fxe <- punif(vx,min,max)}
  if(dist=="weib") {Fxh <- pweibull(vnx,shape,scale);Fxe <- pweibull(vx,shape,scale)}


  lines(vnx,Fxh) 		# theoretical cumulative plot CDF
  legend("bottomright", legend=c("Data","Hyp"), pch=c(1,-1), lty=c(-1,1), merge=T)
  plot(vx,Fx-Fxe,xlab="x", ylab="Diff Empir - Theor"); abline(h=0)
  ret.val <- list(Fx.Emp=Fx,Fx.Theo=Fxh,Fx.Diff=Fx-Fxe)
}

# chisq GOF test with respect to standard normal
# user gives number of equiprobable classes as argument nclass

chisq.gof.norm <- function(x, nclass, param.est){ 
  prob <- rep(1/nclass, nclass)	# equal prob of each class
  expec <- length(x) * prob	# expected number in each class
  numb <- floor(1 + nclass * pnorm(x))	# prob from normal scaled to nclass 
  obser <- tabulate(numb, nclass)	# count observed in each class
  X2 <- sum((expec - obser)^2/expec)	# calc X2
  df <- nclass - 1 - param.est		# df, subtract num par estimated
  p.value <- 1 - pchisq(X2, df)		# pvalue from chisq
  # return result in a list
  ret.val <- list(X2=X2, df=df, p.value=p.value, observed=obser)  
}


invent.mxn <- function(m,n=5,d=1,p,f2="random") {
x <- matrix(nrow=n,ncol=m)
for(i in 1:m){
 if(f2=="random") x[,i] <- rnorm(n,p[i,1],p[i,2]) 
 if(f2=="step"){
  x[,i] <-  seq(p[i,1],p[i,2],(p[i,2]-p[i,1])/(n-1))
 }
}
return(round(x,d))
}


