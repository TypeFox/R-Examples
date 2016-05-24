
#@numb - number of covariates for bb
#@covs - number of total parameters to estimate
#@all the rest in documentation

#Samp implements importance sampling

.samp <- function(t,x,n, Zb, Zw, par, varcv, nsims, keep, numb, covs,
                  erho, esigma, ebeta, ealphab, ealphaw, Rfun){
  import1 <- NULL

  varcv2 <- solve(varcv)/4

  draw <- rmvnorm(nsims, par[covs], varcv2)
  varcv3 <- solve(varcv2)
  phiv <- dmvnorm(draw, par[covs], varcv2, log=T)
  zbmiss <- ifelse(covs[6] == FALSE,TRUE,FALSE)
  zwmiss <- ifelse(covs[(6+numb)] == FALSE, TRUE, FALSE)
  if(zbmiss == TRUE & zwmiss == FALSE){
    draw <- cbind(draw[,1:5], rep(1,nsims), draw[,(5+numb):sum(covs)])
  }
  if(zbmiss == FALSE & zwmiss == TRUE){
    draw <- cbind(draw, rep(1,nsims))
  }
  if(zbmiss == TRUE & zwmiss == TRUE){
    draw <- cbind(draw, rep(1,nsims), rep(1,nsims))
  }
  #for(i in 1:nsims){
      #import1[i] <- -like(as.vector(draw[i,]),
      #t, x, n, Zb, Zw, numb=numb, erho, esigma, ebeta,
      #ealphab, ealphaw, Rfun) - phiv[i]
	#}
  
  #Calculates importance ratio
  import1 <- apply(as.matrix(1:nsims),1,function(i)
                  -like(as.vector(draw[i,]),t, x, n, Zb, Zw,
                  numb=numb, erho, esigma, ebeta, ealphab, ealphaw,Rfun)
                   - phiv[i])

  ok <- !is.nan(import1)
  lnir <- import1-max(import1[ok])
  ir <- NA
  ir[ok] <- exp(lnir[ok])
  #print(mean(is.finite(ir)))
  tst <- ifelse(is.finite(ir), ir>runif(1,0,1), FALSE)
  #print(sum(tst))
  keep <- rbind(keep, draw[tst,])
  return(keep)
}

                                   #@sub -- indeces to create R for
#@Rfun - R function to use
#@bb,bw,sb,sw,rho -- current MLE estimates
#numb numw -- number of covaraites for bb and bw

.createR <- function(sub, Rfun, bb, bw, sb,sw, rho, x, numb, numw){
  out <- NULL 
  lower = cbind(-bb[sub]/sb, -bw[sub]/sw)
  upper = cbind(-bb[sub]/sb+1/sb, -bw[sub]/sw+1/sw)
  mean=c(0,0)
  corr=matrix(c(1,rho,rho,1),nrow=2)

  if (Rfun==1){
    out <- NULL
    makeR <- function (i){
      qi <- pmvnorm(lower=lower[i,], upper=upper[i,], mean=mean,
                    corr=corr)
    }
    out <- foreach(i = 1:length(x[sub]), .combine="c") %dopar% makeR(i)
    #out <- apply(as.matrix(1:length(x[sub])), 1, makeR)
    out <- ifelse(out<0|out==0, 1*10^-322,out)
    out <- log(out)
    #if(sum(is.na(out))>0|sum((out==Inf))>0) print("R not real")
    out <- ifelse(is.na(out)|abs(out==Inf), 999, out)
    return(out)
  }
  if (Rfun==2){
    makeR <- function(i){
      qi <- sadmvn(lower=lower[i,], upper=upper[i,], mean=mean,
                   varcov=corr)
    }
    #out <- foreach(i = 1:length(x[sub]), .combine="c") %dopar% makeR(i)
    out <- apply(as.matrix(1:length(x[sub])), 1, makeR)
    out <- ifelse(out<0|out==0, 1*10^-322,out)
    out <- log(out)
    #if(sum(is.na(out))>0|sum((out==Inf))>0) print("R not real")
    out <- ifelse(is.na(out)|abs(out==Inf), 999, out)
        #return(out)
        # for(i in 1:length(x[sub])){
        # qi <- sadmvn(lower=lower[i,], upper=upper[i,], mean=mean, varcov=corr)
        #qi <- ifelse(qi<0|qi==0, 1*10^-322,qi)
        #out[i] <- log(qi)
        #if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
        #out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
    #}
    return(out)
  }
  if (Rfun==3){
    fun <- function(x) dmvnorm(x,mean,corr)
    for (i in 1:length(x[sub])){
      qi <- adaptIntegrate(fun, lowerLimit=lower[i,],
                           upperLimit=upper[i,])$integral
      out[i] <- log(qi)
      #if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
      out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
    }
    return(out)
  }
  ## THIS VERSION RELIES ON mvnprd FROM PACKAGE mvtnormpcs WHICH IS NOW ARCHIVED
    if (Rfun==4){
        print("Option Rfun==4 is no longer possible")
        stop()
  ##  for(i in 1:length(x[sub])){
  ##    qi <- mvnprd(A=upper[i,], B=lower[i,],
  ##                 BPD=c(rho,rho),INF=rep(2,2))$PROB
  ##    qi <- ifelse(qi<0|qi==0, 1*10^-322,qi)
  ##    out[i] <- log(qi)
      #if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
  ##    out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
  ##  }
  ##  return(out)
  }

if (Rfun==5){
  lower = lower[1,]
  upper = upper[1,]
  #qi <- pmvnorm(lower=lower, upper=upper, mean=mean, corr=corr)
  qi <- sadmvn(lower=lower, upper=upper, mean=mean, varcov=corr)
  qi <- ifelse(qi<1*10^-14, 1*10^-14,qi)
  qi <- log(qi)
  #if(is.na(qi)|abs(qi)==Inf) print ("R not real")
  qi <- ifelse((is.na(qi)|abs(qi)==Inf),999,qi)
  out <- rep(qi,length(x[sub]))
  return(out)
}
}



#@par - solution to likelihood maximization
#@numb - number of covariates for bb
#@covs - number of total parameters to estimate
#@all the rest in documentation

#Samp implements importance sampling

.samp <- function(t,x,n, Zb, Zw, par, varcv, nsims, keep, numb, covs,
                  erho, esigma, ebeta, ealphab, ealphaw, Rfun){
  import1 <- NULL
  varcv2 <- solve(varcv)/4
  draw <- rmvnorm(nsims, par[covs], varcv2)
  varcv3 <- solve(varcv2)
  phiv <- dmvnorm(draw, par[covs], varcv2, log=T)
  #print("samp")
  zbmiss <- ifelse(covs[6] == FALSE,TRUE,FALSE)
  zwmiss <- ifelse(covs[(6+numb)] == FALSE, TRUE, FALSE)
  if(zbmiss == TRUE & zwmiss == FALSE){
    draw <- cbind(draw[,1:5], rep(1,nsims), draw[,(5+numb):sum(covs)])
  }
  if(zbmiss == FALSE & zwmiss == TRUE){
    draw <- cbind(draw, rep(1,nsims))
  }
  if(zbmiss == TRUE & zwmiss == TRUE){
    draw <- cbind(draw, rep(1,nsims), rep(1,nsims))
  }
  #for(i in 1:nsims){
      #import1[i] <- -like(as.vector(draw[i,]),
      #t, x, n, Zb, Zw, numb=numb, erho, esigma, ebeta,
      #ealphab, ealphaw, Rfun) - phiv[i]
	#}
  
  #Calculates importance ratio
  import1 <- apply(as.matrix(1:nsims),1,function(i)
                  -like(as.vector(draw[i,]),t, x, n, Zb, Zw,
                  numb=numb, erho, esigma, ebeta, ealphab, ealphaw,Rfun)
                   - phiv[i])
  ok <- !is.nan(import1)
  lnir <- import1-max(import1[ok])
  ir <- NA
  ir[ok] <- exp(lnir[ok])
  #print(mean(is.finite(ir)))
  tst <- ifelse(is.finite(ir), ir>runif(1,0,1), FALSE)
  #print(sum(tst))
  keep <- rbind(keep, draw[tst,])
  return(keep)
}

#@All arguments are values on the untruncated scale
#This function reparameterizes to the truncated scale

.repar <- function(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw){
  sb=exp(sb0)
  sw=exp(sw0)
  bb=Bb0*(.25+sb^2) + .5 +
    as.matrix(apply(Zb,2, function(x) x - mean(x)))%*%as.matrix(Bb0v)
  bw=Bw0*(.25+sw^2) + .5 +
    as.matrix(apply(Zw,2, function(x) x - mean(x)))%*%as.matrix(Bw0v)
  rho=(exp(2*rho0)-1)/(exp(2*rho0) +1)
  return(c(t(bb), t(bw), sb, sw, rho))
}

#Functions for plotting

#Tomography plot

.tomog <- function(ei.object,title="Tomography Plot with the Data",lci=T){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  #Take out the bounds
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  #Number of precincts n
  n <- dim(bounds)[1]
  #Figure out length of line for LOCOI
  length<-NA
  for(i in 1:n){
    length[i]<-sqrt(abs(bbounds[i,1]-bbounds[i,2])^2 + abs(wbounds[i,1]-wbounds[i,2])^2)
  }
  scale<-((length-min(length))/(max(length)-min(length)))*100
  #Plot
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1),
       col="white", ylab="betaW", xlab="betaB", xaxs="i",
       yaxs="i", main=title)
  if(lci==T){for(i in 1:n){lines(bbounds[i,], wbounds[i,], col=hcl(h=30,c=100,l=scale[i],alpha=1))}}
  else{for(i in 1:n){lines(bbounds[i,], wbounds[i,], col="yellow")}}
}

.tomogd <- function(x,t,n,title,lci=T){
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  n <- dim(bounds)[1]
  #Figure out length of line for LOCOI
  length<-NA
  for(i in 1:n){
      length[i]<-sqrt(abs(bbounds[i,1]-bbounds[i,2])^2 + abs(wbounds[i,1]-wbounds[i,2])^2)
  }
  scale<-((length-min(length))/(max(length)-min(length)))*100
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1),
       col="white", ylab="betaW", xlab="betaB", xaxs="i",
       yaxs="i", main=title)
  if(lci==T){for(i in 1:n){lines(bbounds[i,], wbounds[i,], col=hcl(h=30,c=100,l=scale[i],alpha=1))}}
  else{for(i in 1:n){lines(bbounds[i,], wbounds[i,], col="yellow")}}
}

#Tomography plot with ML contours
.tomogl <- function(ei.object, lci=TRUE){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  Zb <- ei.object$Zb
  Zw <- ei.object$Zw
  phi <- ei.object$phi
  .tomogd(x,t,n, "Tomography Plot with ML Contours", lci=lci)
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]
  Bb0 <- phi[1]
  Bw0 <- phi[2]
  sb0 <- phi[3]
  sw0 <- phi[4]
  rho0 <- phi[5]
  Bb0v <- phi[6:(5+numb)]
  Bw0v <- phi[(6+numb):length(phi)]
  vars <- .repar(Bb0,Bw0,sb0,sw0,rho0, Bb0v, Bw0v, Zb, Zw)
  bb <- vars[1:length(x)]
  bw <- vars[(length(x)+1):(2*length(x))]
  sb <- vars[2*length(x)+1]
  sw <- vars[2*length(x)+2]
  rho <- vars[2*length(x)+3]
  .tomog3 <- function(bb,bw,sb,sw,rho){
    lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),
        centre=c(mean(bb),mean(bw)),level=.914), col="blue",lwd=4)
    lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),
        centre=c(mean(bb),mean(bw)),level=.35), col="red",lwd=4)
    points(mean(bb),mean(bw),col="pink",  pch=15)
  }
  
#tomog4 <- function(bb,bw,sb,sw,rho){
#	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), #scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.914), #col="blue",lwd=1, lty=3)
#	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), #scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.35), #col="red",lwd=1, lty=3)
#	}
#or
#for (i in 1:length(x)){
	#points(mean(bb[i]), mean(bw[i]), col="red", #pch=19, cex=.1)
#	tomog3(bb[i], bw[i], sb, sw, rho)
	#}
  .tomog3(bb,bw,sb,sw,rho)
}

.tomog3 <- function(bb,bw,sb,sw,rho){
  lines(ellipse(matrix(c(1,rho,rho,1),nrow=2),
        scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.914),
        col="blue",lwd=4)
  lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),
        centre=c(mean(bb),mean(bw)),level=.35), col="red",lwd=4)
  points(mean(bb),mean(bw),col="pink",  pch=15)
}


#Tomography plot with 80% CIs
.tomog80CI <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  #Only consider precincts that are heterogeneous
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  .tomogd(x,t,n,"Tomography Plot with 80% CIs",lci=F)
  #Create confidence intervales
  betabcd <- apply(betabs,1,function(x) quantile(x, probs=c(.1,.9)))
  betawcd <- apply(betaws,1,function (x) quantile(x,probs=c(.1,.9)))
  n <- dim(betabcd)[2]
  for(i in 1:n){
    lines(betabcd[,i], sort(betawcd[,i],decreasing=T), col="red",
          lwd=3)
  }
  }
}

#Tomography plot with 95% CIs
.tomog95CI <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  .tomogd(x,t,n,"Tomography Plot with 95% CIs",lci=F)
  betabcd <- apply(betabs,1,function(x) quantile(x,
                                                 probs=c(.025,.975)))
  betawcd <- apply(betaws,1,function (x)
                   quantile(x,probs=c(.025,.975)))
  n <- dim(betabcd)[2]
  for(i in 1:n){
    lines(betabcd[,i], sort(betawcd[,i], decreasing=T), col="red",
          lwd=3)
  }
}
}

#TomogE -- Tomography plot with mean posterior betabs and betaws
.tomogE <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  .tomogd(x,t,n,"Tomography Plot with Mean Posterior Betabs and
Betaws")
  betabm <- apply(betabs,1,mean)
  betawm <- apply(betaws,1,mean)
  points(betabm, betawm, col="red", pch=19)
}
}

#TomogP -- Tomography plot with contours based on mean posterior psi
#tomogd(data$x,data$t, data$n, "Tomography with simulated contours based on mean posterior")
#vars <- apply(psi,2,mean)
#bb <- vars[1]
#bw <- vars[2]
#sb <- vars[3]
#sw <- vars[4]
#rho <- vars[5]
#tomog(bb,bw,sb,sw,rho)
#dev.off()


.tomogP2 <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  psi <- ei.object$psi
  .tomogd(x,t,n,"Tomography Plot with Contours Based on Posterior")
  bbp <- psi[,1:length(x)]
  bwp <- psi[,(length(x)+1):(2*length(x))]
  sbp <- psi[,2*length(x)+1]
  swp <- psi[,2*length(x)+2]
  rhop <- psi[,2*length(x)+3]
  points(mean(bbp), mean(bwp), col="red", pch=19)
  .tomog3(bbp,bwp,mean(sbp),mean(swp),mean(rhop))
#or
#points(mean(bbp), mean(bwp), col="red", pch=19)
#for (i in 1:length(x)){
#	points(mean(bbp[,i]), mean(bwp[,i]), col="red", #pch=19, cex=.1)
#	tomog4(bbp[,i], bwp[,i], mean(sbp), mean(swp), mean #(rhop))
#	}
}
}


#Density plots
#tomogd(data$x,data$t,data$n, "Tomography with contours #from posterior")
#bivn.kde <- kde2d(psi[,1], psi[,2], n=100)
#contour(bivn.kde, add=T, nlevels=7, drawlabels=F)
#points(mean(psi[,1]), mean(psi[,2]), col="red", pch=15)

#Density of betab
.betabd <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)
  betabs <- ei.object$betabs[ok,]
  betabm <- apply(betabs,1,mean)
  plot(density(betabm), xlim=c(0,1),  col="green", xlab="betaB",
       ylab="density across precincts, f(betaB)", main="Density of
betaB")
  vb <- as.vector(betabm)
  for (i in 1:length(vb)){
    lines(c(vb[i], vb[i]), c(0,.25))
  }
}
}

#Density of betaw
.betawd <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betaw)
  betaws <- ei.object$betaws[ok,]
  betawm <- apply(betaws,1,mean)
  plot(density(betawm), xlim=c(0,1), col="green", xlab="betaW",
       ylab="density across precincts, f(betaW)", main="Density of
betaW")
  vw <- as.vector(betawm)
  for (i in 1:length(vw)){
    lines(c(vw[i], vw[i]), c(0,.25))
  }
}
}


#XT plot
.xt <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T
Scatterplot", ylab="T", xlab="X", pch=20)
}

#XTc plot

.xtc <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  circ <- .04
  plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T
Scatterplot with Population Density Circles", ylab="T", xlab="X",
       pch=20)
  minn <- min(n)
  maxn <- max(n)
  for (i in 1:length(x)){
    radius = (n[i]-minn+1)/(1+maxn-minn)
    draw.circle(x[i], t[i], radius*circ)
  }
}
#xtc(data$x,data$t,data$n,.04, "X and T Scatterplot")


#XTfit plot

.xtfit <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  low <- .1
  up <- .9
  circ <- .04
  plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T
Scatterplot with E(T|X) and 80% CIs", ylab="T", xlab="X", pch=20)
  minn <- min(n)
  maxn <- max(n)
  for (i in 1:length(x)){
    radius = (n[i]-minn+1)/(1+maxn-minn)
    draw.circle(x[i], t[i], radius*circ)
  }
  x <- seq(0,1,by=.01)
  betabs <- as.vector(betabs)
  betaws <- as.vector(betaws)
  t <- matrix(ncol=length(x), nrow=length(betabs))
  for(i in 1:length(x)){
    t[,i] <- betabs*x[i] + betaws*(1-x[i])
  }
  et <- apply(t,2,mean)
  lines(x,et, col="yellow")
  lwr <- apply(t,2,function (x) quantile(x, probs=c(low)))
  upr <- apply(t,2,function (x) quantile(x, probs=c(up)))
  lines(x, lwr, col="red")
  lines(x, upr, col="red")
}
}

#xtfit(x,t,n,.04,"",ei.1$betabs, ei.1$betaws, .2,.8)


#XTfitg plot
.xtfitg <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  low <- .1
  up <- .9
  circ <- .04
  plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T
Scatterplot with E(T|X), 80% CIs, and Goodman", ylab="T", xlab="X",
       pch=20)
  minn <- min(n)
  maxn <- max(n)
  for (i in 1:length(x)){
    radius = (n[i]-minn+1)/(1+maxn-minn)
    draw.circle(x[i], t[i], radius*circ)
  }
  x <- seq(0,1,by=.01)
  betabs <- as.vector(betabs)
  betaws <- as.vector(betaws)
  t <- matrix(ncol=length(x), nrow=length(betabs))
  for(i in 1:length(x)){
    t[,i] <- betabs*x[i] + betaws*(1-x[i])
  }
  et <- apply(t,2,mean)
  lines(x,et, col="yellow")
  lwr <- apply(t,2,function (x) quantile(x, probs=c(low)))
  upr <- apply(t,2,function (x) quantile(x, probs=c(up)))
  lines(x, lwr, col="red")
  lines(x, upr, col="red")
  t <- ei.object$t
  x <-ei.object$x
  lm.fit <- lm(t ~ x)
  abline(lm.fit, col="green")
}
}

#Goodman plot
.goodman <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T
Scatterplot with Goodman", ylab="T", xlab="X", pch=20)
  lm.fit <- lm(t ~ x)
  abline(lm.fit, col="red")
}


#Estsims plot

.estsims <- function(ei.object){
  if(!("betabs"%in% names(ei.object))){
   message("Error: This plot function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  colors = runif(length(betabs),26,51)
  plot(betabs, betaws, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i",
       main="Simulations of betaW and betaB", ylab="betaW
simulations", xlab="betaB simulations", pch=20, col=colors, lty=2,
       cex=.25)
}
}

#boundXB


.boundXb <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  truebb <- ei.object$truth[,1]
  bounds <- bounds1(x, t, n)
  plot(x, truebb, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i",
       main="Aggregation Bias for betaB", ylab="True betab", xlab="X",
       pch=20)
  for (i in 1:length(x)){
    lines(c(x[i], x[i]), c(bounds[,1][i], bounds[,2][i]))
  }
  lm.xb <- lm(truebb ~ x)
  abline(lm.xb, lty=2)
}

.boundXw <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  truebw <- ei.object$truth[,2]
  bounds <- bounds1(x, t, n)
  plot(x, truebw, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i",
       main="Aggregation Bias for betaW", ylab="True betaw", xlab="X",
       pch=20)
  for (i in 1:length(x)){
    lines(c(x[i], x[i]), c(bounds[,3][i], bounds[,4][i]))
  }
  lm.xw <- lm(truebw ~ x)
  abline(lm.xw, lty=2)
}


#truth
.truthfn <- function(ei.object){
  n <- ei.object$n
  x <- ei.object$x
  omx <- 1-x
  truebb <- ei.object$truth[,1]
  truebw <- ei.object$truth[,2]
  betabs <- ei.object$betabs
  betaws <- ei.object$betaws
  betab <- ei.object$betab
  betaw <- ei.object$betaw
  truthbb <- sum(truebb*n)/sum(n)
  truthbw <- sum(truebw*n)/sum(n)
  circ=.04
  par(mfrow=c(2,2))
  ag <- .aggs(ei.object)
  plot(density(ag[,1]),
       xlim=c(0,1),ylim=c(0,max(density(ag[,1])$y)+1),
       yaxs="i",xaxs="i", main="Density of Bb Posterior & Truth",
       xlab="Bb",ylab="Density")
  lines(c(truthbb, truthbb), c(0,.25*(max(density(ag[,1])$y)+1)),
        lwd=3)
  plot(density(ag[,2]),
       xlim=c(0,1),ylim=c(0,max(density(ag[,2])$y)+1), yaxs="i",
       xaxs="i", main="Density of Bw Posterior & Truth",
       xlab="Bw",ylab="Density")
  lines(c(truthbw, truthbw), c(0,.25*(max(density(ag[,2])$y)+1)),
        lwd=3)
  plot(betab, truebb, xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",
       xlab="Estimated betab",cex=.1,ylab="True betab")
  minn <- min(x*n)
  maxn <- max(x*n)
  for (i in 1:length(betab)){
    radius = (n[i]*x[i]-minn+1)/(1+maxn-minn)
    draw.circle(betab[i], truebb[i], radius*circ)
  }
  ci80b = .CI80b(ei.object)
  low = mean(abs(ci80b[,1]-betab))
  high = mean(abs(ci80b[,2]-betab))
  abline(0,1)
  lines(c(0,1),c(-low,1-low),lty=2)
  lines(c(0,1),c(high,1+high),lty=2)
  plot(betaw, truebw,
       xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",xlab="Estimated
betaw", ylab="True betaw",cex=.1)
  minn <- min(omx*n)
  maxn <- max(omx*n)
  for (i in 1:length(betaw)){
    radius = (omx[i]*n[i]-minn+1)/(1+maxn-minn)
    draw.circle(betaw[i], truebw[i], radius*circ)
  }
  ci80w = .CI80w(ei.object)
  low = mean(abs(ci80w[,1]-betaw))
  high = mean(abs(ci80w[,2]-betaw))
  abline(0,1)
  lines(c(0,1),c(-low,1-low),lty=2)
  lines(c(0,1),c(high,1+high),lty=2)
}


#Functions for Quantities of Interest

.betaB <- function(ei.object){
    if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ei.object$betab
}
}

.betaW <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ei.object$betaw
}
}

.sbetab <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ei.object$sbetab
}
 }

.sbetaw <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ei.object$sbetaw
}
}
   

.phi <- function(ei.object){
  ei.object$phi
}

.psisims <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ei.object$psi
}
}

.bounds <- function(ei.object){
  bounds1(ei.object$x, ei.object$t, ei.object$n)
}

.abounds <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  omx <- 1-x
  Nb <- n*x
  Nw <- n*omx
  LAbetaB <- weighted.mean(bounds[,1],Nb)
  UAbetaB <- weighted.mean(bounds[,2], Nb)
  LAbetaW <- weighted.mean(bounds[,3], Nw)
  UAbetaW <- weighted.mean(bounds[,4], Nw)
  return(matrix(c(LAbetaB, UAbetaB, LAbetaW,UAbetaW), nrow=2))
}
	
.aggs <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betab <- ei.object$betabs[ok,]
  betaw <- ei.object$betaws[ok,]
  omx <- 1-x
  Nb <- n*x
  Nw <- n*omx
  Bbgg <- vector(mode="numeric", length=dim(betab)[2])
  for (i in 1:dim(betab)[2]){
    Bbgg[i] <- weighted.mean(betab[,i], Nb)
  }
  Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
  for (i in 1:dim(betaw)[2]){
    Bwgg[i] <- weighted.mean(betaw[,i], Nw)
  }
  return(cbind(Bbgg, Bwgg))
}
}
	
.maggs <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betab <- ei.object$betabs[ok,]
  betaw <- ei.object$betaws[ok,]
  omx <- 1-x
  Nb <- n*x
  Nw <- n*omx
  Bbgg <- vector(mode="numeric", length=dim(betab)[2])
  for (i in 1:dim(betab)[2]){
    Bbgg[i] <- weighted.mean(betab[,i], Nb)
  }
  Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
  for (i in 1:dim(betaw)[2]){
    Bwgg[i] <- weighted.mean(betaw[,i], Nw)
  }
  return(c(mean(Bbgg), mean(Bwgg), sd(Bbgg), sd(Bwgg)))
}
}
	
.VCaggs <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betab <- ei.object$betabs[ok,]
  betaw <- ei.object$betaws[ok,]
  omx <- 1-x
  Nb <- n*x
  Nw <- n*omx
  Bbgg <- vector(mode="numeric", length=dim(betab)[2])
  for (i in 1:dim(betab)[2]){
    Bbgg[i] <- weighted.mean(betab[,i], Nb)
  }
  Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
  for (i in 1:dim(betaw)[2]){
    Bwgg[i] <- weighted.mean(betaw[,i], Nw)
  }
  vc <- matrix(c(var(Bbgg), cov(Bbgg, Bwgg), cov(Bbgg, Bwgg),
                 var(Bwgg, Bwgg)), nrow=2)
  return(vc)
}
}
	
.CI80b <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  betab <- ei.object$betabs[ok,]
  lwr <- vector(mode="numeric",length=length(ei.object$x))
  upr <- vector(mode="numeric",length=length(ei.object$x))
  lwr[ok] <- apply(betab, 1, function(x) quantile(x, probs=c(.1)))
  lwr[!ok] <- NA
  upr[ok] <-apply(betab, 1, function(x) quantile(x, probs=c(.9)))
  upr[!ok] <- NA
  return(cbind(lwr,upr))
}
}
	
.CI80w <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  betaw <- ei.object$betaws[ok,]
  lwr <- vector(mode="numeric",length=length(ei.object$x))
  upr <- vector(mode="numeric",length=length(ei.object$x))
  lwr[ok] <- apply(betaw, 1, function(x) quantile(x, probs=c(.1)))
  lwr[!ok] <- NA
  upr[ok] <-apply(betaw, 1, function(x) quantile(x, probs=c(.9)))
  upr[!ok] <- NA
  return(cbind(lwr,upr))
}
}
	
.eaggbias <- function(ei.object){
   if(!("betabs"%in% names(ei.object))){
   message("Error: This eiread function requires an ei.sim object.")
  }
  if("betabs"%in% names(ei.object)){
  x <- ei.object$x
  mbetab <- ei.object$betab
  mbetaw <- ei.object$betaw
  lm.b <- lm(mbetab ~ x)
  lm.w <- lm(mbetaw ~ x)
  output <-
    list(summary(lm.b)$coefficients,summary(lm.w)$coefficients)
  names(output) <- c("betaB", "betaW")
  return(output)
}
}
	
.goodman <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  lm.g <- lm(t ~ x)
  w <- 1-x
  lm.w <- lm(t ~ w)
  BetaW <- summary(lm.g)$coefficients[1,]
  BetaB <- summary(lm.w)$coefficients[1,]
  coefs <- rbind(BetaB, BetaW)
  return(coefs)
}

.movieD <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  n <- dim(bounds)[1]
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white",
       ylab="betaW", xlab="betaB", xaxs="i",yaxs="i",
       main="Tomography Plot")
  input<-0
  while(input!="s"){
    input<-tomogonce(input,last.input)
    last.input<-input
    input<-getinput()
  }
}
getinput <- function(){
  readline("Hit <enter> for next observation, enter observation number, or hit <s> to stop: ")
}

.tomogonce <- function(input,last.input){
  par(mfrow=c(1,1))
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white",
       ylab="betaW", xlab="betaB", xaxs="i",yaxs="i",
       main="Tomography Plot")
  if(input==""){ #input is <enter>, plot next observation
    last.input<-as.integer(last.input)
    input<-last.input+1
    for(i in 1:n){
      lines(bbounds[i,], wbounds[i, ], col="yellow")
    }  
    lines(bbounds[input,], wbounds[input,], col="black")
  }
  else{ #input is observation number
    input<-as.integer(input)
    for(i in 1:n){
      lines(bbounds[i,], wbounds[i,], col="yellow")
    }
    lines(bbounds[input,], wbounds[input,], col="black")
  }
  return(input)
}





#Movie plots

.movieD <- function(ei.object){
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  n <- dim(bounds)[1]
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white",
       ylab="betaW", xlab="betaB", xaxs="i",yaxs="i",
       main="Tomography Plot")
  input<-0
  while(input!="s"){
    input<-.tomogonce(input,last.input, ei.object)
    last.input<-input
    input<-.getinput()
  }
}
.getinput <- function(){
  readline("Hit <enter> for next observation, enter observation number, or hit <s> to stop: ")
}

.tomogonce <- function(input,last.input, ei.object){
    x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  n <- dim(bounds)[1]
  par(mfrow=c(1,1))
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white",
       ylab="betaW", xlab="betaB", xaxs="i",yaxs="i",
       main="Tomography Plot")
  if(input==""){ #input is <enter>, plot next observation
    last.input<-as.integer(last.input)
    input<-last.input+1
    for(i in 1:n){
      lines(bbounds[i,], wbounds[i, ], col="yellow")
    }  
    lines(bbounds[input,], wbounds[input,], col="black")
  }
  else{ #input is observation number
    input<-as.integer(input)
    for(i in 1:n){
      lines(bbounds[i,], wbounds[i,], col="yellow")
    }
    lines(bbounds[input,], wbounds[input,], col="black")
  }
  return(input)
}


.getinput <- function(){
  readline("Hit <enter> for next observation, enter observation number, or hit <s> to stop: ")
}

.postonce<-function(input,last.input,betab,betaw,betabs,betaws){

  if(input==""){ #input is <enter>, will plot next observation
    last.input<-as.integer(last.input)
    input<-last.input+1
  }
  else{ #input is observation number
    input<-as.integer(input) 
  }

  par(mfrow=c(2,2), oma=c(0,0,2,0))

  # plot posterior distribution of precinct parameters
  plot(density(betab[input,]),
    xlim=c(0,1),ylim=c(0,max(density(betab[input,])$y)+1),
    yaxs="i",xaxs="i", main="Posterior Distribution of betaB",
    xlab="Bb",ylab="Density")
  lines(c(0,.25*(max(density(betab[input,])$y)+1)),lwd=3)
  plot(density(betaw[input,]),
    xlim=c(0,1),ylim=c(0,max(density(betaw[input,])$y)+1),
    yaxs="i",xaxs="i", main="Posterior Distribution of betaW",
    xlab="Bw",ylab="Density")
  lines(c(0,.25*(max(density(betaw[input,])$y)+1)),lwd=3)
  
# plot simulated values of betaB and betaW of precinct
  colors = runif(length(betabs),26,51)
  plot(betabs[input,], betaws[input,], xlim=c(0,1), ylim=c(0,1),xaxs="i",
    yaxs="i",main="Simulations of betaW and betaB",
    ylab="betaW simulations", xlab="betaB simulations",
    pch=20,col=colors,lty=2,cex=.25)

  mtext(sprintf("Plots for Observation %d", input),line=0.5,outer=TRUE)

return(input)

}

.movie <- function(ei.object){
  ok <- !is.na(ei.object$betab)
  betabs <- ei.object$betabs[ok,]
  ok <- !is.na(ei.object$betaw)
  betaws <- ei.object$betaws[ok,]

  betab <- ei.object$betabs
  betaw <- ei.object$betaws
  
 input<-1 #initialize at precinct 1
 last.input<-0
 while(input!="s"){
    input<-.postonce(input,last.input,betab,betaw,betabs,betaws)
    last.input<-input
    input<-.getinput()
 }
   
}



#EI RxC Plots

#########
##Heat plot for the eiRxC case
#########

#form = formula
#total = totals for each precinct
#data = data
#names = Names of groups in order: X-axis, Y-axis, Other
#covariate = extra group or covariate to create colors in the plot for
tomogRxC <- function(formula, data, total=NULL, refine=100){
	#require(sp)  Now specified as Import
    noinfocount = 0
	form <- formula
    ##total <- dbuf$total
    #data <- dbuf$data
    dvname <- terms.formula(formula)[[2]]
    covariate <- NA
	#Make the bounds
	rows <- c(all.names(form)[6:(length(all.names(form)))])
    names=rows
    cols <- c(all.names(form)[3])
#    if(sum(data[,rows][,1]<1.1)==length(data[,rows][,1])){
#       data <- round(data*data[,total])
#}
	#print(data[,cols])
    options(warn=-1)
	bnds <- bounds(form, data=data,rows=rows, column =cols,threshold=0)
    options(warn=0)
	#Totals
	dv <- data[, all.names(form)[3]]
	#Assign other category
	bndsoth <- bnds$bounds[[3]]
	oth <- data[,all.names(form)[length(all.names(form))]]

	#Assign x-axis category
	bndsx <- bnds$bounds[[1]]
	xcat <- data[,all.names(form)[6]]

	#Assign y-axis category
	bndsy <- bnds$bounds[[2]]
	ycat <- data[,all.names(form)[7]]

	#Minimums & Maximums
	minx <- bndsx[,1]
	miny <- bndsy[,1]
	minoth <- bndsoth[,1]
	maxx <- bndsx[,2]
	maxy <- bndsy[,2]
	maxoth <- bndsoth[,2]

	#####
	#Starting point when x is at minimum
	##
	#Holding x at its minimum, what are the bounds on y?
	
	#When x is at its minimum, the new dv and total are:
	newdv <- dv - (minx*xcat)
	newtot <- oth + ycat 
	t <- newdv/newtot
	y <- ycat/newtot
	
	#The new bounds on the y category are:

	lby <- cbind(miny, (t - maxoth*oth/newtot)/(y))
	lby[,2] <- ifelse(y==0, 0, lby[,2])
	lowy <- apply(lby,1,max)
	hby <- cbind((t-minoth*oth/newtot)/y,maxy)
	highy <- apply(hby,1,min)

	#####
	#Starting point when x is at maximum
	##
	#Holding x at its maximum, what are the bounds on y?

	#The new bounds on x are:
	newtot <- oth + xcat
	newdv <- dv - (miny*ycat)
	x <- xcat/newtot
	t <- newdv/newtot
	lbx <- cbind(minx, (t-maxoth*oth/newtot)/x)
	lbx[,2] <- ifelse(x==0, 0, lbx[,2])
	lowx <- apply(lbx,1,max)
	hbx <- cbind((t-minoth*oth/newtot)/x,maxx)
	highx <- apply(hbx,1,min)

	#Graph starting points
	#High starting points
	hstr <- cbind(minx, highy)
	#High ending points
	hend <- cbind(highx, miny)

	#Low starting points
	lstr <- cbind(minx, lowy)
	lend <- cbind(lowx, miny)

	xl <- paste("Percent", names[1], dvname[2])
	yl <- paste("Percent", names[2], dvname[2])
	mn <- paste("Tomography Plot in a 2x3 Table (", names[3], " Other Category)*", sep="")

plot(c(0,0), xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i",xlab=xl, ylab=yl, col="white", main=mn)


ok <- !is.na(hstr[,2]) & !is.na(hend[,1])
exp1 <- hstr[ok,2]>=maxy[ok]
exp2 <- hend[ok,1]>=maxx[ok]
exp3 <- lstr[ok,2]<=miny[ok]
exp4 <- lend[ok,1]<=minx[ok]
hstr <- hstr[ok,]
hend <- hend[ok,]
lstr <- lstr[ok,]
lend <- lend[ok,]
dv <- dv[ok]
ycat <- ycat[ok]
oth <- oth[ok]
minoth <- minoth[ok]
xcat <- xcat[ok]
maxy <- maxy[ok]
maxx <- maxx[ok]
contourx <- seq(0,1,by=1/refine)
contoury <- seq(0,1,by=1/refine)
contourz <- matrix(0,nrow=length(contourx), ncol=length(contoury))
for(i in 1:dim(hstr)[1]){
	if((exp1[i] + exp2[i] + exp3[i] + exp4[i])==0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1])
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2])
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        if(area==1){
            noinfocount <- noinfocount+1
        }
		for(j in 1:length(as.vector(contourx))){
			contourz[j,] <- contourz[j,] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs)==1,1+1/(.1+area),0)
        
	}
	}
	if((exp1[i]==1) & (exp2[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(cut)
		kink1y <- c(maxy[i])
		}
	if((exp2[i]==1) & (exp1[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		kink1x <- c(maxx[i])
		kink1y <- c(cut)
	}
	if((exp2[i]==1 & exp1[i]==1)){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(maxx[i], cut2)
		kink1y <- c(cut, maxy[i])
	}
	if((exp3[i]==1) & (exp4[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(cut)
		kink2y <- c(miny[i])
		}
	if((exp4[i]==1) & (exp3[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		kink2x <- c(minx[i])
		kink2y <- c(cut)
	}
	if((exp3[i]==1 & exp4[i]==1)){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(minx[i], cut2)
		kink2y <- c(cut, miny[i])
	}
	if((exp3[i] + exp4[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
	    c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        if(area==1){
            noinfocount <- noinfocount+1
        }
		for(j in 1:length(as.vector(contourx))){
			contourz[j,] <- contourz[j,] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs)==1,1+1/(.1+area),0)
         
	}
		
	}
	if((exp1[i] + exp2[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1])
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2])
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        if(area==1){
            noinfocount <- noinfocount+1
        }
		for(j in 1:length(as.vector(contourx))){
			contourz[j,] <- contourz[j,] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs)==1,1+1/(.1+area),0)
	}	}
	if((exp1[i] + exp2[i])!=0 & (exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        if(area==1){
            noinfocount <- noinfocount+1
        }
		for(j in 1:length(as.vector(contourx))){
			contourz[j,] <- contourz[j,] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs)==1,1+1/(.1+area),0)
	}
	}
	}
image(contourx[2:refine], contoury[2:refine], contourz[2:refine,2:refine], , col=sort(heat.colors(refine), decreasing=T), xlab=xl, ylab=yl, main=mn, xlim=c(0,1), add=T)
for(i in 1:dim(hstr)[1]){
	if((exp1[i] + exp2[i] + exp3[i] + exp4[i])==0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1])
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2])
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
	}
	if((exp1[i]==1) & (exp2[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(cut)
		kink1y <- c(maxy[i])
		}
	if((exp2[i]==1) & (exp1[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		kink1x <- c(maxx[i])
		kink1y <- c(cut)
	}
	if((exp2[i]==1 & exp1[i]==1)){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(maxx[i], cut2)
		kink1y <- c(cut, maxy[i])
	}
	if((exp3[i]==1) & (exp4[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(cut)
		kink2y <- c(miny[i])
		}
	if((exp4[i]==1) & (exp3[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		kink2x <- c(minx[i])
		kink2y <- c(cut)
	}
	if((exp3[i]==1 & exp4[i]==1)){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(minx[i], cut2)
		kink2y <- c(cut, miny[i])
	}
	if((exp3[i] + exp4[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
	    c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
	}
	if((exp1[i] + exp2[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1])
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2])
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.25)
	}
	if((exp1[i] + exp2[i])!=0 & (exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
        polygon(xaxs,yaxs, col=rgb(0,0,0,alpha=0), border="black", lty=1, lwd=.225)
	}
	}
message(paste("There are", noinfocount, "tomography polygons with no information"))
#contour(contourx[2:refine], contoury[2:refine], contourz[2:refine,2:refine], nlevels=nrow(data), method="simple", col="black", add=TRUE, vfont = c("sans serif", "plain"), drawlabels=F)
par(xpd=NA)
text(.78,-.18, paste("*There are", noinfocount, "tomography polygons with no information"), cex=.75)
par(xpd=FALSE)
	}


############
##Bounds plot for eiRxC case
############

#form = formula
#total = totals for each precinct
#data = data
#names = Names of groups in order: X-axis, Y-axis, Other
#covariate = extra group or covariate to create colors in the plot for
.bndplot <- function(dbuf){
    #require(sp)  Now specified as Import
	form <- dbuf$formula
    total <- dbuf$total
    data <- dbuf$data
    n <- nrow(data)
    print(n)
    covariate <- NA
	#Make the bounds
	rows <- c(all.names(form)[6:(length(all.names(form)))])
    names=rows
    cols <- c(all.names(form)[3])
    if(sum(data[,rows][,1]<1.1)==length(data[,rows][,1])){
        data <- round(data*data[,total])
}
	bnds <- bounds(form, data=data, rows=rows, column =cols,threshold=0)
    
	#Totals
	dv <- data[, all.names(form)[3]]

	#Assign other category
	bndsoth <- bnds$bounds[[3]]
	oth <- data[,all.names(form)[length(all.names(form))]]

	#Assign x-axis category
	bndsx <- bnds$bounds[[1]]
	xcat <- data[,all.names(form)[6]]

	#Assign y-axis category
	bndsy <- bnds$bounds[[2]]
	ycat <- data[,all.names(form)[7]]

	#Minimums & Maximums
	minx <- bndsx[,1]
	miny <- bndsy[,1]
	minoth <- bndsoth[,1]
	maxx <- bndsx[,2]
	maxy <- bndsy[,2]
	maxoth <- bndsoth[,2]

	#####
	#Starting point when x is at minimum
	##
	#Holding x at its minimum, what are the bounds on y?
	
	#When x is at its minimum, the new dv and total are:
	newdv <- dv - (minx*xcat)
	newtot <- oth + ycat 
	t <- newdv/newtot
	y <- ycat/newtot
	
	#The new bounds on the y category are:
	lby <- cbind(miny, (t - maxoth*oth/newtot)/(y))
	lby[,2] <- ifelse(y==0, 0, lby[,2])
	lowy <- apply(lby,1,max)
	hby <- cbind((t-minoth*oth/newtot)/y,maxy)
	highy <- apply(hby,1,min)

	####
	#Starting point when x is at maximum
	##
	#Holding x at its maximum, what are the bounds on y?
	#The new bounds on x are:
	newtot <- oth + xcat
	newdv <- dv - (miny*ycat)
	x <- xcat/newtot
	t <- newdv/newtot
	lbx <- cbind(minx, (t-maxoth*oth/newtot)/x)
	lbx[,2] <- ifelse(x==0, 0, lbx[,2])
	lowx <- apply(lbx,1,max)
	hbx <- cbind((t-minoth*oth/newtot)/x,maxx)
	highx <- apply(hbx,1,min)

	###
	#Graph starting points
	####

	#High starting points
	hstr <- cbind(minx, highy)
	#High ending points
	hend <- cbind(highx, miny)

	#Low starting points
	lstr <- cbind(minx, lowy)
	lend <- cbind(lowx, miny)
	
	#Colors for covariates
	if(!is.na(covariate)){
		redg <- data[,covariate]/(oth-data[,covariate])/max(data[,covariate]/(oth-data[,covariate]))
		blug <- 1-redg
	}
	if(is.na(covariate)){
		redg <- rep(.5, length(minx))
		blug <- rep(.5, length(minx))	
	}
	
	#Graph labels
	xl <- paste("Percent", names[1], "Vote Democrat")
	yl <- paste("Percent", names[2], "Vote Democrat")
	mn <- paste("Tomography Plot in a 2x3 Table (", names[3], " Other Category)", sep="")
	
	#Initial plot
	plot(c(0,0), xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i",xlab=xl, ylab=yl, col="white", main=mn)

	#Only non-NA starting pts are OK
	ok <- !is.na(hstr[,2]) & !is.na(hend[,1])
	
	#All different types of polygons
	exp1 <- hstr[ok,2]>=maxy[ok]
	exp2 <- hend[ok,1]>=maxx[ok]
	exp3 <- lstr[ok,2]<=miny[ok]
	exp4 <- lend[ok,1]<=minx[ok]
	
	#Subsets all the variables
	hstr <- hstr[ok,]
	hend <- hend[ok,]
	lstr <- lstr[ok,]
	lend <- lend[ok,]
	dv <- dv[ok]
	ycat <- ycat[ok]
	oth <- oth[ok]
	minoth <- minoth[ok]
	xcat <- xcat[ok]
	maxy <- maxy[ok]
	maxx <- maxx[ok]
	redg <- redg[ok]
	blug <- blug[ok]
	
for(i in 1:dim(hstr)[1]){
	#4 corner polygon
	if((exp1[i] + exp2[i] + exp3[i] + exp4[i])==0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1])
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2])
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
		#if(area>0 & !is.nan(area)){alpha = 1/(1+b*area)}
		#if(area>0 & !is.nan(area)){alpha = .5-.5*area}
		#if(area>0 & !is.nan(area)){alpha = ((1-area)^(3)+.2)*.83}
		if(area>0 & !is.nan(area)){alpha = min(.5/(area*(n)),1)}
		if(area==0 | is.nan(area)){alpha=.05}
		border = alpha 
		polygon(xaxs,yaxs, col=rgb(redg[i],0,blug[i],alpha=alpha), border=rgb(redg[i],0,blug[i],alpha=1), lty=2)
	}
	
	#Create more corners
	if((exp1[i]==1) & (exp2[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(cut)
		kink1y <- c(maxy[i])
		}
	if((exp2[i]==1) & (exp1[i])==0){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		kink1x <- c(maxx[i])
		kink1y <- c(cut)
	}
	if((exp2[i]==1 & exp1[i]==1)){
		cut <- (dv[i]-(oth[i])*minoth[i])/ycat[i] - maxx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*minoth[i])/xcat[i] - maxy[i]*ycat[i]/xcat[i]
		kink1x <- c(maxx[i], cut2)
		kink1y <- c(cut, maxy[i])
	}
	if((exp3[i]==1) & (exp4[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(cut)
		kink2y <- c(miny[i])
		}
	if((exp4[i]==1) & (exp3[i])==0){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		kink2x <- c(minx[i])
		kink2y <- c(cut)
	}
	if((exp3[i]==1 & exp4[i]==1)){
		cut <- (dv[i]-(oth[i])*maxoth[i])/ycat[i] - minx[i]*xcat[i]/ycat[i]
		cut2 <- (dv[i]-(oth[i])*maxoth[i])/xcat[i] - miny[i]*ycat[i]/xcat[i]
		kink2x <- c(minx[i], cut2)
		kink2y <- c(cut, miny[i])
	}
	
	#Plot 5-sided polygon
	if((exp3[i] + exp4[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
	c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
		#if(area>0 & !is.nan(area)){alpha = 1/(1+b*area)}
		if(area>0 & !is.nan(area)){alpha = min(.5/(area*(n)),1)}
		#if(area>0 & !is.nan(area)){alpha = ((1-area)^(3)+.2)*.53}
		if(area==0 | is.nan(area)){alpha=.05}
		border = alpha 
		polygon(xaxs,yaxs, col=rgb(redg[i],0,blug[i],alpha=alpha), border=rgb(redg[i],0,blug[i],alpha=1), lty=2)
		
	}
	
	#Another 5-sided polygon
	if((exp1[i] + exp2[i])==0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1])
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2])
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
		if(area>0 & !is.nan(area)){alpha = min(.5/(area*(n)),1)}
		#if(area>0 & !is.nan(area)){alpha = ((1-area)^(3)+.2)*.53}
		if(area==0 | is.nan(area)){alpha=.05}
		border = alpha 
		polygon(xaxs,yaxs, col=rgb(redg[i],0,blug[i],alpha=alpha), border=rgb(redg[i],0,blug[i],alpha=1), lty=2)
	}
	
	#Plot 6-sided polygons
	if((exp1[i] + exp2[i])!=0 & (exp3[i] + exp4[i])!=0){
		xaxs <- c(hstr[i,1],  lstr[i,1],kink2x,lend[i,1],hend[i,1], kink1x)
		xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
		yaxs <- c(hstr[i,2], lstr[i,2],kink2y,lend[i,2], hend[i,2], kink1y)
		yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
		side1 <- c(xaxs,xaxs[1])
		side2 <- c(yaxs,yaxs[1])
		c1 <- sum(side1[1:(length(side1)-1)]*side2[2:length(side2)])
		c2 <- sum(side1[2:(length(side1))]*side2[1:(length(side2)-1)])
		area <- abs(c1-c2)/2
		#if(area>0 & !is.nan(area)){alpha = 1/(1+b*area)}
		if(area>0 & !is.nan(area)){alpha = min(.5/(area*(n)),1)}
		#if(area>0 & !is.nan(area)){alpha = ((1-area)^(3)+.2)*.53}	
		if(area==0 | is.nan(area)){alpha=.05}
		border = alpha 
		polygon(xaxs,yaxs, col=rgb(redg[i],0,blug[i],alpha=alpha), border=rgb(redg[i],0,blug[i],alpha=1), lty=2)
		}
	}
}


	

