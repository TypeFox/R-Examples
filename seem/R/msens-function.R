# function to calculate sens metrics for multiple runs

msens <- function(t.X, param, fileout, pdfout=F, resp=F) {

# sensitivity analysis
# arguments:
# t.X time seq and model results
# param list of plabel, param, and factor values
# fileout prefix to store file and plots

 # number of runs and labels for runs
  npar <- length(param$plab)
  if(npar==1) nruns <- length(param$pval) else  nruns <- dim(param$pval)[1]

 # labels for runs
  if(npar==1) {
    valrun <- param$pval;labrun <- param$plab   
    param$pval <- matrix(param$pval)  
  } else {
    valrun <- array(); labrun <- param$plab[1] 
    for(i in 2:npar) labrun <- paste(labrun,param$plab[i],sep=";")
    for(j in 1:nruns){
     valrun[j] <- param$pval[j,1]
     for(i in 2:npar) valrun[j] <- paste(valrun[j],param$pval[j,i],sep=";")  
     }
  }
  legrun <- paste(labrun,"=",valrun,sep="")

 # labels for metrics and responses
 ndex=6; label.m <- c("Diff Max", "Diff Min", "Diff End", "Avg Diff", "RMS Diff", "Max Diff")
 rdex <- 4; label.r <- c("Max", "Min", "End", "Avg")

 # parameters
 p <- param$pval; pnom <- param$pnom
 # finds position of nominal in the set of parameter values 
 for(i in 1:nruns) if(sum(p[i,]-pnom)==0) nnom=i 
 # % change in parameter (abs is needed for negative values)
 pp <- p; for(i in 1:npar) pp[,i] <- 100*(p[,i] - p[nnom,i])/abs(p[nnom,i])
 # parameters as standard variables
 ps <- pp; for(i in 1:npar) ps[,i] <- 100*(pp[,i] - mean(pp[,i]))/sd(pp[,i])

 # time and variables
 t <- t.X[,1]; X <- t.X[,-1]; nt <- length(t)

 # create arrays and set to zero
 d.max <- array()
 d.min <- array()
 d.end <- array()
 avg.d <- array()
 err.m <- matrix(nrow=nt,ncol=nruns) 
 rms.d <- array()
 max.d <- array()

 # for each parameter value calculate the metrics
 for(j in 1:nruns){
   d.max[j] <- max(X[,j])
   d.min[j] <- min(X[,j])
   d.end[j] <- X[nt,j]
   avg.d[j] <- mean(X[,j])
   err.m[,j] <- X[,j] - X[,nnom]
   rms.d[j] <- sqrt(mean((err.m[,j])^2))
   max.d[j] <- max(abs(err.m[,j]))
}
 y <- cbind(d.max, d.min, d.end, avg.d, rms.d, max.d)

 # relative error and relative metrics
 err.r <- err.m
 for(j in 1:nruns){
   for(i in 1:nt){
    if (X[i,nnom] ==0) err.r[i,j] <- 0 else
      err.r[i,j] <- (err.m[i,j])/X[i,nnom]
   }
   rms.d[j] <- sqrt(mean((err.r[,j])^2))
   max.d[j] <- max(abs(err.r[,j]))
 }
 yp <- cbind(d.max, d.min, d.end, avg.d, rms.d, max.d)

 # % change in metrics
 yn <- y[nnom,]
  for(i in 1:rdex){
 # checks for all zero values for a metric
   if(yn[i] == 0) yp[,i]=0
 # calculate % if nonzero
   else yp[,i] <- 100*(y[,i]-yn[i])/yn[i] 
  }
  for(i in (rdex+1):ndex) yp[,i] <- 100*(yp[,i]) 

 Ma <- data.frame(p,round(y,2))
 names(Ma) <-  c(param$plab,label.m)
 Mp <- data.frame(pp,round(yp,2))
 names(Mp) <-  c(param$plab,label.m)

 if(param$fact == T) {
 
 # pairwise selection for more than 2 param
 nc <- factorial(npar)/(factorial(2)*factorial(npar-2))
 combo <- matrix(nrow=nc,ncol=2);k=1
 for(i in 1:(npar-1)){
  for(j in (i+1):npar) {combo[k,]<- c(i,j);k=k+1}
 }
 # pairwise response matrices 
 nlev <- dim(param$pv)[1]
 x <- list(); xp <- list()
 j=0
 for(k in 1:nc){ 
 for(i in 1:ndex){
  j=j+1 
  x[[j]] <-  round(tapply(y[,i],list(p[,combo[k,1]],p[,combo[k,2]]),mean),3)
  xp[[j]] <- round(tapply(yp[,i],list(p[,combo[k,1]],p[,combo[k,2]]),mean),3)
 }
 }
 
 if(resp==F){ # only for factorial sensitivity
  # ANOVA 
  aov.pvalue <- matrix(nrow=npar,ncol=ndex); AOV.pvalue <- list()
  friedman.pvalue <- matrix(nrow=1,ncol=ndex); Friedman.pvalue <- list()

  for(k in 1:nc){
  for(i in 1:ndex){
   aov.pvalue[,i] <- summary(aov(y[,i]~p[,combo[k,1]]+p[,combo[k,2]],data=Ma))[[1]][1:2,5]
   friedman.pvalue[,i] <- friedman.test(y[,i], p[,combo[k,1]],p[,combo[k,2]])[3]$p.value
  }
  AOV.pvalue[[k]] <- data.frame(round(aov.pvalue,3))
  names(AOV.pvalue[[k]]) <- label.m
  row.names(AOV.pvalue[[k]]) <- param$plab
  Friedman.pvalue[[k]] <- data.frame(round(friedman.pvalue,3))
  names(Friedman.pvalue[[k]]) <- label.m
  #row.names(Friedman.pvalue[[k]]) <- param$plab

 }

 Eff <- matrix(NA,2,ndex); Eff.Int <- Eff 
 for(i in 1:ndex){
 Y <- abs(xp[[i]])
 Eff[1,i] <- mean(Y[-2,])- mean(Y[2,])
 Eff[2,i] <- mean(Y[,-2])- mean(Y[,2])
 Eff.Int[1,i] <-  mean(Y[-2,-2]) - mean(Y[-2,2])
 Eff.Int[2,i] <-  mean(Y[-2,-2]) - mean(Y[2,-2])
 }
 Eff <- data.frame(round(Eff,3)); names(Eff) <- label.m; row.names(Eff)<- param$plab
 Eff.Int <- data.frame(round(Eff.Int,3)); names(Eff.Int) <- label.m; row.names(Eff.Int)<- param$plab
 } # end fact sensitivity no resp

 } else { # only for sampling
 # regression
 ycoeff <- matrix(nrow=ndex,ncol=npar); yp.est <- yp;r2 <- array()
 for(i in 1:rdex){
  ycoeff[i,] <- lm(yp[,i]~0+ps)$coeff
  r2[i] <- summary(lm(yp[,i]~0+ps))$r.square
 }
 for(i in (rdex+1):ndex){
  ppa <- abs(ps)
  ycoeff[i,] <- lm(yp[,i]~0+ppa)$coeff
  r2[i] <- summary(lm(yp[,i]~0+ppa))$r.square
 }
 for(k in 1:nruns){
  for(i in 1:4) yp.est[k,i] <- sum(ycoeff[i,]*ps[k,]) 
  for(i in 5:6) yp.est[k,i] <- sum(ycoeff[i,]*ppa[k,])
 } 
 slope <- t(ycoeff)
 reg.slope.R2 <- data.frame(round(rbind(slope, t(r2)),3))
 names(reg.slope.R2) <- label.m
 row.names(reg.slope.R2) <- c(param$plab,"R2")
 } # end sampling

 # graphics
 if(pdfout==T) pdf(paste(fileout,".pdf",sep=""))
 mat <- matrix(1:6,3,2,byrow=T)
 nf <- layout(mat, widths=rep(7/2,2), heights=rep(7/3,3), TRUE)
 par(mar=c(4,4,1,.5),xaxs="r",yaxs="r")

 if(param$fact == T) {
  if(resp==F){ # only for factorial sensitivity
  xt <- length(param$pv[,2])
  x1 <- p[,1]; x2 <- p[,2]
  for(j in 1:ndex){
   interaction.plot(x1,x2,y[,j],ylab=label.m[j],xlab=param$plab[1],
                 type="b", pch=1:xt, trace.label=param$plab[2])
  }
  x1 <- pp[,1]; x2 <- pp[,2]
  for(j in 1:ndex){
   interaction.plot(x1,x2,yp[,j],ylab=label.m[j],xlab=param$plab[1],
                 type="b", pch=1:xt, trace.label=param$plab[2])
  }
 }else {# only for factorial response resp
   mat <- matrix(1:4,2,2,byrow=T)
   nf <- layout(mat, widths=rep(7/3,2), heights=rep(7/3,2), TRUE)
   par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")

   label.r.nc <-  rep(label.r,nc)
   for(k in 1:nc){
   x1 <- abs(param$pv[,combo[k,1]]); x2 <- abs(param$pv[,combo[k,2]])
   for(j in (rdex*(k-1)+1):(rdex*k)) {
    image(x1,x2,x[[j]],xlab=param$plab[combo[k,1]],ylab=param$plab[combo[k,2]],col=grey(10:20/20),pty="s")
    contour(x1,x2,x[[j]],add=T)
    title(label.r.nc[j])
   } # j loop
  } # k loop
 } # else resp

 } else { # only for sampling
  for(j in 1:ndex){
   plot(yp[,j],yp.est[,j],ylab=label.m[j],xlab="Metric %")
   abline(a=0,b=1)
  }
 } 
 if(pdfout==T) dev.off()
 
  if(param$fact == T & resp==F) return(list(Ma=Ma,Mp=Mp, AOV.pvalue=AOV.pvalue, Friedman.pvalue=Friedman.pvalue, Eff=Eff, Eff.Int=Eff.Int))
  if(param$fact == T & resp==T) return(list(Ma=Ma,Mp=Mp,x=x))
  if(param$fact == F) return(list(Ma=Ma,Mp=Mp, reg.slope.R2=reg.slope.R2))

}

