## File IDPSurvival/R/isurvdiff.r
##
## IDPSurvival package for R (http://www.R-project.org)
##
## Copyright (C) 2014 Dalle Molle Institute for Artificial Intelligence.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

isurvdiff <- function(formula, data, groups=c(1,2), s=0.25,
                    alternative = c("two.sided", "less", "greater"),
                    exact=NULL,  level = 0.95, display=TRUE, 
                    nsamples=10000, rope=0, tmax=NULL) {
  Call <- match.call()
  Call[[1]] <- as.name('survfit')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
  indx <- match(c('formula', 'data', 'weights'), names(Call), nomatch=0)
  #It's very hard to get the next error message other than malice
  #  eg survfit(wt=Surv(time, status) ~1) 
  if (indx[1]==0) stop("a formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  m <- eval.parent(temp)
  
  Terms <- terms(formula)
  ord <- attr(Terms, 'order')
  if (length(ord) & any(ord !=1))
    stop("Interaction terms are not valid for this function")
  
  if (!length(groups)==2)
    stop("Can only handle 2 groups")
  
  n <- nrow(m)
  Y <- model.extract(m, 'response')
  if (!is.Surv(Y)) stop("Response must be a survival object")
  

  ll <- attr(Terms, 'term.labels')
  if (length(ll) == 0) {X <- factor(rep(1,n))  # ~1 on the right
  } else {
    X <- strata(m[ll])
    groups <- paste(ll,groups,sep='=')
  }
  
  if (!is.Surv(Y)) stop("y must be a Surv object")
  
  
  if (attr(Y, 'type') == "right") {
    temp <- do.call(".isurvdiff.test",
                    list(X,Y,groups,s,alternative,exact,level,display,nsamples,rope,tmax))
  } else {
    stop("Can only handle right censored data")
  }
  return(temp)
}


.isurvdiff.test <- function(x, y, groups=c(1,2), s=0.25,
                    alternative = c("two.sided", "less", "greater"),
                    exact,  level = 0.95, display=TRUE, 
                    nsamples=10000, rope=0, tmax=NULL) {
  
  alternative <- match.arg(alternative)
  
  if (!is.Surv(y)) stop("y must be a Surv object")
  if (!is.factor(x)) stop("x must be a factor")
  if (attr(y, 'type') != 'right')
    stop("Can only handle right censored data")
  
  # Each of the necessary output objects is originally a list with one
  #  element per strata.  
  time <- vector('list', 2)
  survUP <- vector('list', 2)
  n.cens <- vector('list', 2)
  n.event<- vector('list', 2)
  n.risk <- vector('list', 2)
  Pxy <- vector('list', 2)
  g.opt <- vector('list', 2)
  strata <- integer(2)
  
  uniquex <- sort(unique(x))
  check <- match(groups,uniquex)
  
  if (!is.na(match( NA,check)))
    stop("incorrect group names")

  casewt=rep(1,length(x))
  for (i in 1:2) {
    who <- (x== groups[i])
    
    temp <- tapply(casewt[who], 
                   list(factor(y[who,1]), 
                        factor(y[who,2], levels=0:1)), sum)
    temp <- ifelse(is.na(temp), 0, temp)
    ntemp  <- (dim(temp))[1]

    time[[i]] <- type.convert(dimnames(temp)[[1]], as.is=TRUE,
                              dec=getOption("OutDec"))
    nevent <- as.vector(temp[,2])
    ncens  <- as.vector(temp[,1])
    nrisk  <- rev(cumsum(rev(nevent+ncens)))
    terms <- (s+nrisk-nevent)/(s+nrisk)
    terms[which(nrisk==0)]<-1 #avoid 0/0 cases
   
    survUP[[i]] <- cumprod(terms)
    strata[i] <- ntemp    
    n.event[[i]] <- nevent
    n.cens[[i]]  <- ncens
    n.risk[[i]]  <- nrisk
  }
  
  if (missing(exact) || is.null(exact)) {
    if (strata[1]<100|strata[2]<100) { exact<-TRUE }
    else { 
      exact<-FALSE
      nsamples<-200000
    }
  } else if (exact==FALSE) {nsamples<-200000}
  
# TEST
  if(is.null(tmax)) D.sUP <- list(-diff(c(1,survUP[[1]])),-diff(c(1,survUP[[2]])))
  else D.sUP <- list(-diff(c(1,survUP[[1]],0)),-diff(c(1,survUP[[2]],0))) 

  time.mat <- list(matrix(time[[1]],1,strata[1]),matrix(time[[2]],1,strata[2]))
  ix.ref <- c(1,2)
  iy.ref <- c(2,1)

  objfun <- function(g, s, dsup, nevent, ncens, nrisk, ntemp, A) {
    gful <- rep(0,ntemp)
    gful[c(1,which(ncens>0)+1)] <- g
    G <- s-cumsum(gful)+gful
    pterm <- (nrisk-nevent+G-gful)/(nrisk+G)
    pterm[which((nrisk+G)==0)] <- 1
    pterm[ntemp]<-0
    dslow <- -diff(c(1,cumprod(pterm)))
#     A<-matrix(rep(1,(ntemp)*(ntemp)),ntemp,ntemp)
#     A[lower.tri(A, diag = FALSE)]<-0
#     diag(A) <- 1/2
    EP <- sum((dsup%*%A)*dslow)
    EP
  } 

  confun <- function(g, s, dsup, nevent, ncens, nrisk, ntemp,A) { #nrisk=nriskm
    con <- sum(g)
    con
  } 

  for (i in 1:2) { #i=1: upper P(G1<G2) i.e. lower P(G2<G1), quindi prendo upper G2 e lower G1; 
                   #i=2: lower P(G1<G2) quindi prendo upper G1 e lower G2; 
                   # qui chiamo x e y in modo che risulti lower P(y<x)
    
    ix <- ix.ref[i]
    iy <- iy.ref[i]
    ntemp <- strata[ix]
    ncens <- n.cens[[ix]]
    dsup <- D.sUP[[iy]]
    x <-  time.mat[[ix]]
    y <- t(time.mat[[iy]])
    
    epsilon <- min(c(diff(y),diff(x)))/1000
    ng <- length(ncens[ncens>0])+1
    ind <- 1:ntemp+cumsum(c(0,ncens[-ntemp]>=1))+1
    
    # add 'prior data'
    ncens.app <- rep(0,ng+ntemp)
    ncens.app[ind]<-ncens
    nevent.app <- rep(0,ng+ntemp)
    nevent.app[ind]<-n.event[[i]]
    nrisk.app <- rep(0,ng+ntemp)
    nrisk.app[ind]<-n.risk[[i]]
    nrisk.app[which(ncens.app>0)+1] <- nrisk.app[which(ncens.app>0)]-ncens.app[which(ncens.app>0)]
    nrisk.app[1] <- nrisk.app[2]
    x.app <- rep(0,ng+ntemp)
    x.app[ind]<-x
    x.app[which(ncens.app>0)+1] <- x.app[which(ncens.app>0)]+epsilon
    x.app <- matrix(x.app,1,ntemp+ng)
   
    n.x <- length(x.app)
    n.y <- length(y)
    Xmat <- .repmat(x.app,n.y,1)
    Ymat <- .repmat(y,1,n.x)
    A <- (sign(Xmat-Ymat)+1)/2 

	if(!is.null(tmax)) {
      A[(Ymat>tmax)*(Xmat>tmax)] <- 0.5
      A <-rbind(A,0.5*(x.app>tmax))     
    }

    g0 <- rep(s/ng,ng)
    if (s > 0) {
      g <- solnp(g0, objfun, eqfun=confun, eqB=s,  LB = rep(0,ng), UB = rep(s,ng), control=list(trace=0),
                    s=s, ncens=ncens.app, nevent=nevent.app, dsup=dsup, nrisk=nrisk.app, ntemp=ntemp+ng, A=A)
      g.opt[[i]] <- g$pars
    } else {
      g <- list("pars"=0*g0)                                                                                              
      g.opt[[i]] <- 0*g0
    }

    if(exact==TRUE)
    { 
      #sample Sx
      t.g <- c(1,which(ncens.app>0)+1)
      t.cens <- which(ncens.app>0)
      t.event <- which(nevent.app>0)
      a=c(nevent.app[t.event], ncens.app[t.cens],g$pars)
      Hk <- rdirichlet(nsamples,a)
      Hevent <- .repmat(matrix(0,1,ng+ntemp),nsamples,1)
      Hevent[,t.event] <- Hk[,1:length(t.event)] 
      Hcens <- .repmat(matrix(0,1,ng+ntemp),nsamples,1)
      Hcens[,t.cens] <- Hk[,length(t.event)+1:length(t.cens)] 
      Hevent[,t.g] <- Hevent[,t.g]+Hk[,length(t.event)+length(t.cens)+1:length(t.g)] 
      Hrisk <- 1-t(apply(Hevent+Hcens, 1, cumsum))+Hevent+Hcens
      Hrisk[which(Hrisk<0)]=0
      terms <- (Hrisk-Hevent)/Hrisk
#       terms[which(Hrisk==0)]<-1 # avoids 0/0 
      terms[,ng+ntemp]<-0 # the last element in terms is 0/Hrisk[,ng+ntemp]
                          #(Hrisk[,ng+ntemp] may go to 0 but is never =0)      
      SsampLOW <- cbind(matrix(1,nsamples,1), t(apply(terms, 1, cumprod)))
      dslow <- -t(apply(SsampLOW, 1, diff))
  
      #Sample Sy
      y.cens <- n.cens[[iy]] 
      y.event<- n.event[[iy]]
      t.cens <- which(y.cens>0)
      t.event <- which(y.event>0)
      a=c(y.event[t.event], y.cens[t.cens], s)
      Hk <- rdirichlet(nsamples,a)
      Hevent <- .repmat(matrix(0,1,n.y+1),nsamples,1)
      Hevent[,t.event] <- Hk[,1:length(t.event)] 
      Hevent[,n.y+1] <- Hevent[,n.y+1]+Hk[,length(a)]
      Hcens <- .repmat(matrix(0,1,n.y+1),nsamples,1)
      Hcens[,t.cens] <- Hk[,length(t.event)+1:length(t.cens)]
      Hrisk <- 1-t(apply(Hevent+Hcens, 1, cumsum))+Hevent+Hcens
      Hrisk[which(Hrisk<0)]=0
      terms <- (Hrisk-Hevent)/Hrisk
      terms[which(Hrisk==0)]<-1 # avoids 0/0

	  if(is.null(tmax)) SsampUP <- cbind(matrix(1,nsamples,1), t(apply(terms[,-(n.y+1)], 1, cumprod)))
	  else SsampUP <- cbind(matrix(1,nsamples,1), t(apply(terms[,-(n.y+1)], 1, cumprod)),matrix(0,nsamples,1)) 

      dsup <- -t(apply(SsampUP, 1, diff))
      
      #Compute P(y<x)
      dataXY <- (dsup%*%A)*dslow
      Pxy[[i]] <- apply(dataXY,1,sum)
    }  
    else
    {
      #Use Normal approximation
      
      #compute mean
      meanP <- objfun(g$par, s=s, ncens=ncens.app, nevent=nevent.app, 
                      dsup=dsup, nrisk=nrisk.app, ntemp=ntemp+ng, A=A)
      
      #compute variance
      y.risk <- c(n.risk[[iy]]+s,s)
      y.event<- c(n.event[[iy]],s)
      y.event[which(y.risk==0)] <- 1                                                                                        
      y.risk[which(y.risk==0)] <- 1                                                                                         
      n.y<-n.y+1
      Wi <- c(1, cumprod((y.risk-y.event+1)/(y.risk+1)))
      Wi <-  y.event/(y.risk+1)*Wi[-(n.y+1)]
      Wi <- matrix( Wi,n.y,1)
      Wj <- c(1, cumprod((y.risk-y.event)/(y.risk)))
      Wj <-  y.event/(y.risk)*Wj[-(n.y+1)]
      Wj <- matrix(Wj,1,n.y)
      WW <- Wi%*%Wj
      Wii <- c(1,cumprod((y.risk-y.event+1)*(y.risk-y.event)/((y.risk+1)*(y.risk))))
      Wii <-  y.event*(1+y.event)/((y.risk+1)*(y.risk))*Wii[-(n.y+1)]
      Wii <- matrix(Wii,n.y,1)
      WW[lower.tri(WW, diag = FALSE)] <- t(WW)[lower.tri(WW, diag = FALSE)]
      diag(WW) <- Wii
      
      t.g <- c(1,which(ncens.app>0)+1)
      g.app <- rep(0,n.x)
      g.app[t.g] <- g$pars
      x.event <- nevent.app + g.app
      x.risk <- nrisk.app + rev(cumsum(rev(g.app)))
      x.event[which(x.risk==0)] <- 1                                                                                        
      x.risk[which(x.risk==0)] <- 1                                                                                         
      Vi <- c(1, cumprod((x.risk-x.event+1)/(x.risk+1)))
      Vi <- x.event/(x.risk+1)*Vi[-(n.x+1)]
      Vi <- matrix(Vi,n.x,1)
      Vj <- c(1, cumprod((x.risk-x.event)/(x.risk)))
      Vj <- x.event/x.risk*Vj[-(n.x+1)]
      Vj <- matrix(Vj,1,n.x)
      VV <- Vi%*%Vj
      Vii <- c(1,cumprod((x.risk-x.event+1)*(x.risk-x.event)/((x.risk+1)*(x.risk))))
      Vii <- x.event*(1+x.event)/((x.risk+1)*x.risk)*Vii[-(n.x+1)]
      Vii <- matrix(Vii,n.x,1)
      VV[lower.tri(VV, diag = FALSE)] <- t(VV)[lower.tri(VV, diag = FALSE)]
      diag(VV) <- Vii
      
      if(is.null(tmax)) A <- t(cbind(t(A),rep(0,n.x)))
      VarP <- sum(diag(t(A)%*%WW%*%A%*%VV))-meanP^2;
      
      Pxy[[i]]  <- rnorm(nsamples,meanP,sqrt(VarP))
      
    }
  }

  data12l <- Pxy[[2]]
  data12u <- 1-Pxy[[1]]

  alphav<-1-level

  #compute credible interval for the lower
  SD <- sort(data12l);
  xll = SD[ceiling(alphav*nsamples/2)] #left bound
  xlr = SD[ceiling((1-alphav/2)*nsamples)] #right bound
  
  #compute credible interval for the upper
  SD = sort(data12u);
  xul = SD[ceiling(alphav*nsamples/2)] #left bound
  xur = SD[ceiling((1-alphav/2)*nsamples)] #right bound
  
  cred_Lower_bound<-c(xll,xlr)
  cred_Upper_bound<-c(xul,xur)
  
  
  
  switch(alternative,
         "two.sided" = {
           if((xlr<0.5-rope && xur<0.5-rope) || (xll>0.5+rope && xul>0.5+rope)){
             h <- 1; #hypotheis H1 is accepted
           }else{
             if((xll<0.5-rope && xlr>=0.5+rope) && (xul<0.5-rope && xur>=0.5+rope)) 
               h <- 0 #hypotheis H0 is accepted
             else
               h <- 2 #no enough information to discrimnate hypotheses H0  and H1
           }
           prob <-NULL
         },
         "greater" = {
           areal <- length(which(data12l>0.5+rope))/nsamples;
           areau <- length(which(data12u>0.5+rope))/nsamples;
           #              
           if (areal>1-alphav && areau>1-alphav){
             h<-1 #hypotheis H1 is accepted
           } else {
             if((areal<=1-alphav && areau>1-alphav) ||  (areal>1-alphav && areau<=1-alphav))
               h<-2 #no enough information to discriminate hypotheses H0  and H1
             else
               h<-0  #hypotheis H0 is accepted
           }
           prob<-c(areal,areau)
         },
         "less" = {
           areal <- length(which(1-data12l>0.5+rope))/nsamples;
           areau <- length(which(1-data12u>0.5+rope))/nsamples;
           #              
           if (areal>1-alphav && areau>1-alphav){
             h<-1 #hypotheis H1 is accepted
           } else {
             if((areal<=1-alphav && areau>1-alphav) ||  (areal>1-alphav && areau<=1-alphav))
               h<-2 #no enough information to discriminate hypotheses H0  and H1
             else
               h<-0  #hypotheis H0 is accepted
           }
           prob<-c(areau,areal)
         })
  
  if (display==TRUE)  {
#     plot.new()    
    ## calculate the density - don't plot yet
    densLower <- density(data12l)
    densUpper <- density(data12u)
    ## calculate the range of the graph
    xlim <- c(0,1)
    ylim <- range(0,densUpper$y, densLower$y)
    #pick the colours
    LowerCol <- rgb(1,0,0,0.2)
    UpperCol <- rgb(0,0,1,0.2)
    LowInt<- rgb(1,0,0,0.8)
    UpInt <- rgb(0,0,1,0.8)
    
    
    if (alternative=="greater")
    {    ## plot the Lowers and set up most of the plot parameters
      plot(densLower, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
           panel.first = grid(),col=LowInt, axes=FALSE)
      xlabel  <- seq(0, 1, by = 0.1)
      axis(1, at = xlabel, las = 1)
      axis(2)
      box()
      lines(densUpper, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
            panel.first = grid(),col=UpInt)
      #put our density plots in
      dens<- densLower
      polygon(c(0.5, dens$x[dens$x>.5 & dens$x <= 1], xlim[2]), c(0, dens$y[dens$x>=.5 & dens$x <= 1], 0),density = -1, col = LowerCol)
      dens<- densUpper
      polygon(c(0.5, dens$x[dens$x>.5 & dens$x <= 1], xlim[2]), c(0, dens$y[dens$x>=.5 & dens$x <= 1], 0),density = -1, col = UpperCol)
      ## add a legend in the corner
      legend('topleft',c('Lower Distribution','Upper Distribution'),
             fill = c(LowerCol, UpperCol), bty = 'n',
             border = NA)
      #       get( getOption( "device" ) )()
      usr <- par( "usr" )
      height <- usr[ 4 ]-usr[ 3 ] 
      width <- usr[ 2 ]-usr[ 1 ]
      text( usr[ 1 ]+width/50, usr[ 4 ]-height/5, paste("Area Lower=",areal),     adj = c( 0, 1 ), col = LowInt )
      text( usr[ 1 ]+width/50,usr[ 4 ]-1.25*height/5, paste("Area Upper=",areau),    adj = c( 0, 1 ), col = UpInt )
      
    }
    if (alternative=="less")
    {    ## plot the Lowers and set up most of the plot parameters
      plot(densLower, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
           panel.first = grid(),col=LowInt, axes=FALSE)
      xlabel  <- seq(0, 1, by = 0.1)
      axis(1, at = xlabel, las = 1)
      axis(2)
      box()
      lines(densUpper, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
            panel.first = grid(),col=UpInt)
      #put our density plots in
      dens<- densLower
      polygon(c(xlim[1], dens$x[dens$x>.0 & dens$x <= 0.5], 0.5), c(0, dens$y[dens$x>=.0 & dens$x <= 0.5], 0),density = -1, col = LowerCol)
      dens<- densUpper
      polygon(c(xlim[1], dens$x[dens$x>.0 & dens$x <= 0.5], 0.5), c(0, dens$y[dens$x>=.0 & dens$x <= 0.5], 0),density = -1, col = UpperCol)
      ## add a legend in the corner
      legend('topleft',c('Lower Distribution','Upper Distribution'),
             fill = c(LowerCol, UpperCol), bty = 'n',
             border = NA)
      #       get( getOption( "device" ) )()
      usr <- par( "usr" )
      height <- usr[ 4 ]-usr[ 3 ] 
      width <- usr[ 2 ]-usr[ 1 ]
      text( usr[ 1 ]+width/50, usr[ 4 ]-height/5, paste("Area Lower=",areal),     adj = c( 0, 1 ), col = LowInt )
      text( usr[ 1 ]+width/50, usr[ 4 ]-1.25*height/5, paste("Area Upper=",areau),    adj = c( 0, 1 ), col = UpInt )
    }
    
    if (alternative=="two.sided")
    {    ## plot the Lowers and set up most of the plot parameters
      plot(densLower, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
           panel.first = grid(),col=LowInt, axes=FALSE)
      xlabel  <- seq(0, 1, by = 0.1)
      axis(1, at = xlabel, las = 1)
      axis(2)
      box()
      lines(densUpper, xlim = xlim, ylim = ylim, xlab = 'P(X<Y)', ylab='Probability', main='IDP rank-sum test',
            panel.first = grid(),col=UpInt)
      #put our density plots in
      dens<- densLower
      polygon(c(xll, dens$x[dens$x>xll & dens$x <= xlr], xlr), c(0, dens$y[dens$x>=xll & dens$x <= xlr], 0),density = -1, col = LowerCol)
      dens<- densUpper
      polygon(c(xul, dens$x[dens$x>xul & dens$x <= xur], xur), c(0, dens$y[dens$x>=xul & dens$x <= xur], 0),density = -1, col = UpperCol)
      ## add a legend in the corner
      legend('topleft',c('Lower Distribution','Upper Distribution'),
             fill = c(LowerCol, UpperCol), bty = 'n',
             border = NA)
      #       get( getOption( "device" ) )()
      usr <- par( "usr" )
      height <- usr[ 4 ]-usr[ 3 ] 
      width <- usr[ 2 ]-usr[ 1 ]
      xll<- round(xll*100)/100
      xul<- round(xul*100)/100
      xlr<- round(xlr*100)/100
      xur<- round(xur*100)/100
      text( usr[ 1 ]+width/50, usr[ 4 ]-height/5, paste("Cred Lower Int=[",xll,",",xlr,"]"),     adj = c( 0, 1 ), col = LowInt )
      text( usr[ 1 ]+width/50, usr[ 4 ]-1.25*height/5, paste("Cred Upper Int=[",xul,",",xur,"]"),    adj = c( 0, 1 ), col = UpInt )
    }
   
    
  }
  names(strata)<-groups
  est<-list("h"=h,"prob"=prob,"Lower.Cred.Int"=cred_Lower_bound,"Upper.Cred.Int"=cred_Upper_bound,"alternative"=alternative,"strata"=strata,"exact"=exact)
  class(est) <- "isurvdiff"
  return(est)
}

print.isurvdiff <- function(x, ...)
{
  x$Lower.Cred.Int[1]<- round(x$Lower.Cred.Int[1]*1000)/1000
  x$Lower.Cred.Int[2]<- round(x$Lower.Cred.Int[2]*1000)/1000
  x$Upper.Cred.Int[1]<- round(x$Upper.Cred.Int[1]*1000)/1000
  x$Upper.Cred.Int[2]<- round(x$Upper.Cred.Int[2]*1000)/1000
  h<-x$h
  cat("\n\n---------------------------------------------\n")
  cat("Result of the IDP RANK-SUM hypothesis test\n")
  cat(paste("h =",h," --> ")) 
#   cat("\n")
  if (x$h==2)
    cat("Indeterminate\n")
  if (!is.null(x$prob)){
    if (x$alternative=="greater"){
      if (x$h==1)
        cat("Y is greater than X\n")
      if (x$h==0)
        cat("Y is not greater than X\n")
      cat("\nLower Probability of the hypothesis: ")
      cat(x$prob[1])
      cat("\nUpper Probability of the hypothesis: ")
      cat(x$prob[2])}
    if (x$alternative=="less"){
      if (x$h==1)
        cat("Y is lower than X\n")
      if (x$h==0)
        cat("Y is not lower than X\n")
      cat("\nLower Probability of the hypothesis: ")
      cat(x$prob[1])
      cat("\nUpper Probability of the hypothesis: ")
      cat(x$prob[2])}
      
  } else {
    if (x$h==1)
      cat("Y and X are different\n")
    if (x$h==0)
      cat("Y and X are not different\n")
  }
  cat("\n")
  cat("\nLower Central credible intervals\n")
  cat(paste("[",x$Lower.Cred.Int[1],",",x$Lower.Cred.Int[2],"]"))
  cat("\nUpper Central credible intervals\n")
  cat(paste("[",x$Upper.Cred.Int[1],",",x$Upper.Cred.Int[2],"] \n"))
  cat("---------------------------------------------")
}





