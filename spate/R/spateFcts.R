get.real.dft.mat <- function(wave,indCos,ns=4,n){
  c1 <- sqrt(n*n/2)
  c2 <- sqrt(n*n)
  base <- rep(0,dim(wave)[2])
  phi <- matrix(0,n*n,dim(wave)[2])
  for(n2 in 0:(n-1)){
    for(n1 in 0:(n-1)){
      i <- n2*n+n1+1
      base[1:ns] <- cos(c(n1,n2)%*%(wave[,1:ns]/n))/c2
      base[indCos] <- cos(c(n1,n2)%*%(wave[,indCos]/n))/c1
      base[indCos+1] <- sin(c(n1,n2)%*%(wave[,indCos+1]/n))/c1
      phi[i,] <- base
    }
  }
  return(phi)
}

matern.spec <- function(wave,n,ns=4,rho0,sigma2,nu=1,norm=TRUE){
  if(nu<=0){
    print("Error: nu needs to be positive")
    return()
  }
  if(sigma2<=0){
    print("Error: sigma2 needs to be positive")
    return()
  }
  NF <- dim(wave)[2]
  w2 <- apply(wave^2,2,sum)
  d <- 2
  specnum <- (2^(nu-1))*nu * ((1/rho0)^(2*nu))
  specdenom <- (pi^(d/2))*((1/rho0)^2 + w2)^(nu + d/2)
  spec<- specnum/specdenom
  spec[1:ns] <- spec[1:ns]/2
  if(norm){
    spec <- spec*(n*n)/sum(spec)
  }else{
    spec <- spec/sum(spec)
  }
  spec=sigma2*spec
  return(spec)
}
innov.spec <- function(wave,n,ns=4,rho0,sigma2,zeta,rho1,alpha,gamma,nu=1,dt=1,norm=TRUE){
  if(nu<0){
    print("Error: nu needs to be positive")
    return()
  }
  if(sigma2<0){
    print("Error: sigma2 needs to be positive")
    return()
  }
  if(zeta<0){
    print("Error: zeta needs to be positive")
    return()
  }
  if(gamma<0){
    print("Error: gamma needs to be positive")
    return()
  }
  if(alpha<0 & alpha>=(pi/2)){
    print("Error: alpha needs to be between 0 and pi/2")
    return()
  }
  spec <- matern.spec(wave=wave,n=n,ns=ns,rho0=rho0,sigma2=1,nu=nu,norm=norm)
  NF <- dim(wave)[2]
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }else{
    Tr <- cbind(c(cos(alpha),-gamma*sin(alpha)),c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Tr)%*%Tr)
  }
  DiffDamp <- -apply(wave*Sig%*%wave,2,sum)-rep(zeta,NF)
  spec <- sigma2*spec*(1-exp(2*dt*DiffDamp))/-2/DiffDamp
  return(spec)
}

get.propagator <- function(wave,indCos,zeta,rho1,gamma,alpha,muX,muY,dt=1,ns=4){
  ##if(zeta!=Inf){}
  NF <- dim(wave)[2]
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }else{
    Tr <- cbind(c(cos(alpha),-gamma*sin(alpha)),c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Tr)%*%Tr)
  }
  DiffDamp <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,NF)
  Adv <- dt*c(muX,muY)%*%wave
  G <- matrix(0,ncol=NF,nrow=NF)
  diag(G) <- c(exp(DiffDamp[1:ns]),exp(DiffDamp[(ns+1):NF])*cos(Adv[(ns+1):NF]))
  diag(G[indCos,indCos+1]) <- -exp(DiffDamp[indCos])*sin(Adv[indCos])
  diag(G[indCos+1,indCos]) <- exp(DiffDamp[indCos])*sin(Adv[indCos])
  return(G)
}

propagate.spectral <- function(alphat,spateFT=NULL,n=NULL,Gvec=NULL,par=NULL){
  if(is.null(spateFT)) spateFT <- spate.init(n=n,T=1)
  if(is.null(Gvec)){
    if(is.null(par)){
      print("Either 'Gvec' or 'par' needs to be given.")
      return()
    }else{
      Gvec <- get.propagator.vec(wave=spateFT$wave,indCos=spateFT$indCos,zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],muX=par[7],muY=par[8],dt=1,ns=spateFT$ns)
    }
  }
  xtp1 <- .C("propagate_spectral",xtp1=as.double(rep(0,spateFT$n*spateFT$n)),as.double(alphat),as.double(Gvec$G11C),as.double(Gvec$G11),as.double(Gvec$G12),as.integer(length(spateFT$indCos)),as.integer(spateFT$ns))$xtp1
  return(xtp1)
}
get.propagator.vec <- function (wave, indCos, zeta, rho1, gamma, alpha, muX, muY, dt = 1, ns = 4) {
  NF <- dim(wave)[2]
  if (rho1 == 0) {
    Sig <- cbind(c(0, 0), c(0, 0))
  }
  else {
    Tr <- cbind(c(cos(alpha), -gamma * sin(alpha)), c(sin(alpha), gamma * cos(alpha)))/rho1
    Sig <- solve(t(Tr) %*% Tr)
  }
  DiffDamp <- -dt * apply(wave * Sig %*% wave, 2, sum) - dt * 
    rep(zeta, NF)
  Adv <- dt * c(muX, muY) %*% wave
  G11C <- exp(DiffDamp[1:ns])
  G11 <- exp(DiffDamp[indCos]) * cos(Adv[indCos])
  G12 <- -exp(DiffDamp[indCos]) * sin(Adv[indCos])
  G <- list(G11C=G11C,G11=G11,G12=G12)
  return(G)
}

real.fft <- function(w,n,inv=TRUE,indFFT=NULL){
  if(is.null(indFFT)){
    indFFT <- index.complex.to.real.dft(n)
  }
  wh <- .C("real_fft",n=as.integer(n),wh=as.double(w),inverse=as.integer(sum(inv)),indCos=as.integer(indFFT$indCos),indW=as.integer(indFFT$indW),indWCon=as.integer(indFFT$indWCon), NFc=as.integer(length(indFFT$indCos)))$wh
  return(wh)
}


## real.fftR <- function(y,n,inv=TRUE,indCos=NULL,indW=NULL,indWCon=NULL){
##   if(is.null(indCos) | is.null(indW)){
##     waveC <- waveR <- array(0,c(2,n*n))
##     idx <- c(rep(0:(n/2),n/2+1),rep(1:(n/2-1),n/2-1))
##     idy <- c(as.vector(apply(matrix(0:(n/2)),1,rep,times=(n/2+1))),as.vector(apply(matrix((-n/2+1):(-1)),1,rep,times=(n/2-1))))
##     wa <- rbind(idx,idy)
##     sinex=c(which(idx==0 & idy ==0),which(idx==(n/2) & idy ==0),which(idx==0 & idy ==(n/2)),which(idx==(n/2) & idy ==(n/2)))
##     indCos <- 2*(1:(dim(wa)[2]-4))-1+4
##     if(is.null(indW)){
##       waveR[,1:4] <- wa[,sinex]
##       waveR[,indCos] <- wa[,-sinex]
##       waveR[,indCos+1] <- wa[,-sinex]
##       waveR[2,waveR[2,]<0] <- waveR[2,waveR[2,]<0]+n
##       waveC <- t(cbind(rep(0:(n-1),n),as.vector(apply(matrix(0:(n-1)),1,rep,times=n))))  
##       indW <- match(waveR[1,]+1i*waveR[2,],waveC[1,]+1i*waveC[2,])
##     }
##     if(is.null(indWCon) & !inv){
##       waveCon <- (n-waveR[,indCos])%%n
##       indWCon <- match(waveCon[1,]+1i*waveCon[2,],waveC[1,]+1i*waveC[2,])
##     }
##   }
##   if(inv){
##     ym <- matrix(y,ncol=n)
##     ## ymh <- fft(t(ym),inverse=TRUE)/n
##     ## fftV <- as.vector(t(ymh))
##     ymh <- fft(ym,inverse=TRUE)/n##same thing
##     fftV <- as.vector(ymh)
    
##     fftV <- fftV[indW]
##     yR <- y
##     yR[1:4] <- fftV[1:4]
##     yR[indCos] <- sqrt(2)*Re(fftV[indCos])
##     yR[indCos+1] <- sqrt(2)*Im(fftV[indCos+1])
##     yR <- as.numeric(Re(yR))
##   }else{
##     fftV <- y
##     fftV[indW[1:4]] <- y[1:4]
##     fftV[indW[indCos]] <- (y[indCos]+1i*y[indCos+1])/sqrt(2)
##     fftV[indWCon] <- (y[indCos]-1i*y[indCos+1])/sqrt(2)
##     ##yR <- as.vector(Re(t(fft(t(matrix(fftV,ncol=n)),inverse=FALSE))))/n
##     yR <- as.vector(Re(fft(matrix(fftV,ncol=n),inverse=FALSE)))/n
##   }
##   return(yR)
## }
  
real.fft.TS <- function(w,n,T,inv=TRUE,indFFT=NULL){
  if(class(w)=="matrix"){
    w <- TSmat.to.vect(w)
    Mat <- TRUE
  }else{
    Mat <- FALSE
  }
  if(length(w)!=(n*n*T)){
    print("Error: 'w' needs to be a vector of length n*n*T or a matrix of dimension T x n*n.")
    return()
  }
  if(is.null(indFFT)){
    indFFT <- index.complex.to.real.dft(n)
  }
  wh <- .C("TSreal_fft",n=as.integer(n),T=as.integer(T),wh=as.double(w),inverse=as.integer(sum(inv)), indCos=as.integer(indFFT$indCos), indW=as.integer(indFFT$indW), as.integer(indFFT$indWCon), as.integer(length(indFFT$indCos)))$wh
  if(Mat) wh <- vect.to.TSmat(wh,T=T)
  return(wh)
}

vnorm <- function(v) return(sqrt(sum(v^2)))
cols <- function() return( c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000"))
Prho0 <- function(rho0,log=FALSE){
  if(rho0>100) c <- 0 else c <- 1
  if(log) return(log(c)) else return(c)
}
Psigma2 <- function(sigma2,log=FALSE){
  if(sigma2<0){
    if(log) return(log(0)) else return(0)
  }else{
    if(log) return(-1/2*log(sigma2)) else return(sigma2^(-1/2))
  }
}
Pgamma <- function(gamma,log=FALSE){
  if(gamma<1/100 | gamma > 100){
    if(log) return(log(0)) else return(0)
  }else{
    if(log) return(-log(gamma)) else return(gamma^(-1))
  }
}
Pzeta <- function(zeta,log=FALSE){
    if(zeta>=0) c <- 1 else c <- 0
    if(log) return(log(c)) else return(c)
}
Palpha <- function(alpha,log=FALSE){
  if((0)<=alpha & alpha<=(pi/2)){
    if(log) return(0) else return(1)
  }else{
    if(log) return(-Inf) else return(0)
  }
}
Prho1 <- function(rho1,log=FALSE){
  if(rho1>100) c <- 0 else c <- 1
  if(log) return (log(c)) else return(c)
}
Pmux <- function(mux,log=FALSE){
  if(-0.5<=mux & mux<=0.5){
    if(log) return(0) else return(1)
  }else{
    if(log) return(-Inf) else return(0)
  }
}
Pmuy <- function(muy,log=FALSE){
  if(-0.5<=muy & muy<=0.5){
    if(log) return (0) else return (1)
  }else{
    if(log) return(-Inf) else return(0)
  }
}
Ptau2 <- function(tau2,log=FALSE){
   if(tau2<0){
    if(log) return(log(0)) else return(0)
  }else{
    if(log) return(-1/2*log(tau2)) else return(tau2^(-1/2))
  }
}
Plambda <- function(lambda,log=FALSE){
    if(lambda>=0) c <- 1 else c <- 0
    if(log) return(log(c)) else return(c)
}

wave.numbers <- function(n){
  if(n%%2==1){
    print("Error: n must be even")
    return()
  }else{
    idx <- c(rep(0:(n/2),n/2+1),rep(1:(n/2-1),n/2-1))
    idy <- c(as.vector(apply(matrix(0:(n/2)),1,rep,times=(n/2+1))),as.vector(apply(matrix((-n/2+1):(-1)),1,rep,times=(n/2-1))))
    wa <- rbind(idx,idy)
    sinex=c(which(idx==0 & idy ==0),which(idx==(n/2) & idy ==0),which(idx==0 & idy ==(n/2)),which(idx==(n/2) & idy ==(n/2)))
    indCos <- 2*(1:(dim(wa)[2]-4))-1+4
    wave <- array(0,c(2,n*n))
    wave[,1:4] <- wa[,sinex]
    wave[,indCos] <- wa[,-sinex]
    wave[,indCos+1] <- wa[,-sinex]
    wave <- wave*2*pi##/n
    return(list(wave=wave,indCos=indCos))
  }
}

index.complex.to.real.dft <- function(n){
  waveC <- waveR <- array(0,c(2,n*n))
  idx <- c(rep(0:(n/2),n/2+1),rep(1:(n/2-1),n/2-1))
  idy <- c(as.vector(apply(matrix(0:(n/2)),1,rep,times=(n/2+1))),as.vector(apply(matrix((-n/2+1):(-1)),1,rep,times=(n/2-1))))
  wa <- rbind(idx,idy)
  sinex=c(which(idx==0 & idy ==0),which(idx==(n/2) & idy ==0),which(idx==0 & idy ==(n/2)),which(idx==(n/2) & idy ==(n/2)))
  indCos <- 2*(1:(dim(wa)[2]-4))-1+4
  waveR[,1:4] <- wa[,sinex]
  waveR[,indCos] <- wa[,-sinex]
  waveR[,indCos+1] <- wa[,-sinex]
  waveR[2,waveR[2,]<0] <- waveR[2,waveR[2,]<0]+n
  waveC <- t(cbind(rep(0:(n-1),n),as.vector(apply(matrix(0:(n-1)),1,rep,times=n))))  
  indW <- match(waveR[1,]+1i*waveR[2,],waveC[1,]+1i*waveC[2,])
  waveCon <- (n-waveR[,indCos])%%n
  indWCon <- match(waveCon[1,]+1i*waveCon[2,],waveC[1,]+1i*waveC[2,])
  return(list(indCos=indCos,indW=indW,indWCon=indWCon))
}

TSmat.to.vect <- function(mat){
  return(as.vector(t(mat)))
}
vect.to.TSmat <- function(vect,T=1){
  return(matrix(vect,nrow=T,byrow=TRUE))
}
spate.sim <- function(par,n,T,seed=NULL,StartVal=NULL,nu=1){
  if(length(par)!=9){
    print("Error: 'par' needs to be of length 9")
    return()
  }
  spateFT <- spate.init(n=n,T=T)
  spec <- innov.spec(wave=spateFT$wave,n=n,ns=spateFT$ns,rho0=par[1],sigma2=par[2],zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],nu=nu,dt=1,norm=TRUE)
  Gvec <- get.propagator.vec(wave=spateFT$wave,indCos=spateFT$indCos,zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],muX=par[7],muY=par[8],dt=1,ns=spateFT$ns)
  if(!is.null(seed)) set.seed(seed)
  Innov <- vect.to.TSmat(rep(sqrt(spec),T)*rnorm(n*n*T),T=T)
  if(!is.null(StartVal)){
    Innov[1,] <- real.fft(StartVal,n,inv=TRUE)
  }
  alpha <- Innov
  for(t in 2:T){
    alpha[t,] <- propagate.spectral(alphat=alpha[t-1,],spateFT=spateFT,Gvec=Gvec)+Innov[t,]
  }
  xi <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alpha),n=n,T=T,inv=FALSE),T=T)
  w <- xi+matrix(rnorm(n*n*T,sd=sqrt(par[9])),nrow=T)
  names(par) <- c("rho0","sigma2", "zeta","rho1","gamma","alpha","muX","muY","tau2")
  ret <- list(xi=xi,w=w,alpha=alpha,par=par,n=n,T=T)
  class(ret) <- "spateSim"
  return(ret)
}
plot.spateSim <- function(x,...,plotXi=TRUE,plotW=FALSE){
  if(plotXi) spate.plot(xi=x$xi,...)
  if(plotW) spate.plot(xi=x$w,...)
}
print.spateSim  <- function(x,...){
  print(paste("spateSim object with T=",x$T,", n=",x$n," and parameters:",sep=""))
  print(signif(x$par,digits=3))
}
summary.spateSim <- function(object,...){
  print(paste("spateSim object with T=",object$T,", n=",object$n," and parameters:",sep=""))
  print(signif(object$par,digits=3))
}
spate.plot <- function(xi, nx=NULL, whichT=NULL, format = "ImgTogether", ToFile = FALSE, path=NULL, file=NULL, indScale = FALSE , main=NULL, mfrow = NULL, imagesize = c(1000, 1000), zlim = NULL, breaks = NULL,...){
  if(is.null(whichT)) whichT <- 1:dim(xi)[1]
  if(is.null(zlim)) zlim <- c(min(xi[whichT,]), max(xi[whichT,]))
  if(is.null(breaks)) breaks = seq(from = zlim[1],to = zlim[2], length.out = length(cols())+1)
  if(is.null(nx)){
    nx <- ny <- sqrt(dim(xi)[2])
  }else{
    ny <- dim(xi)[2]/nx
  }
  plot.xi <- function(){
    for(i in whichT){
      if(is.null(main)){
        maint <- paste("t=",i,sep="")
      }else if(length(main)==length(whichT)){
        maint <- main[i]
      }else if(length(main)==1){
        maint <- main
      }else print("'main' must either be NULL or a character vector of length equal to the number of time points or 1.")
      if(format=="ImgSeparate" & ToFile) jpeg(paste(path,file,t,".jpeg",sep=""),width = imagesize[1], height = imagesize[2])
      if(indScale) image(1:nx,1:ny,matrix(xi[i,],nrow=nx),col=cols(),main=maint,...)
      else image(1:nx,1:ny,matrix(xi[i,],nrow=nx),col = cols(),zlim=zlim,breaks=breaks,main=maint,...)
      if(format=="ImgSeparate" & ToFile) dev.off()
    }
  }
  if(format=="ImgTogether"){
    if(ToFile) jpeg(paste(path,file,".jpeg",sep=""),width = imagesize[1], height = imagesize[2])
    if(is.null(mfrow)){
      if(length(whichT)==1) mfrow=c(1,1)
      else if(length(whichT)==2) mfrow=c(1,2)
      else if(length(whichT)<=4) mfrow=c(2,2)
      else if(length(whichT)<=6) mfrow=c(2,3)
      else if(length(whichT)<=9) mfrow=c(3,3)
      else if(length(whichT)<=12) mfrow=c(3,4)
      else if(length(whichT)<=16) mfrow=c(4,4)
      else mfrow=c(5,5)
    }
    par(mfrow=mfrow,...)
    plot.xi()
    if(ToFile) dev.off()
  }else if(format=="ImgSeparate"){
    plot.xi()
  }
}



spate.init <- function(n,T,NF=n*n){
  waveL <- wave.numbers(n)
  waveL$ns <- 4
  IndFour <- 1:(n*n)
  if(NF<(n*n)){
    cut <- quantile(apply(waveL$wave,2,vnorm)[order(apply(waveL$wave,2,vnorm))],probs=NF/n/n)
    IndFour <- which(apply(waveL$wave,2,vnorm)<=cut)
    if(length(IndFour)!=NF) print(paste("Warning: ",length(IndFour)," Fourier functions (instead of ",NF,") are used since when using only ",NF," there is anisotropy in the reduced dimensional basis.",sep="")) 
    waveL$wave <- waveL$wave[,IndFour]
    waveL$ns <- 4-sum(is.na(match(1:4,IndFour)))
    waveL$indCos <- match(IndFour[match(waveL$indCos,IndFour)[!is.na(match(waveL$indCos,IndFour))]],IndFour)
    NF <- length(IndFour)
   }
  indFFT <- index.complex.to.real.dft(n)
  spateFT <- list(n=n,T=T,wave=waveL$wave,indCos=waveL$indCos,ns=waveL$ns,IndFour=IndFour,indFFT=indFFT)
  class(spateFT) <- "spateFT"
  return(spateFT)
}


ffbs <- function(y,lp,G,Sigma,H,Omega,N=dim(y)[2],T=dim(y)[1],NF=dim(G)[1],lglk=FALSE,BwSp=TRUE,filt=FALSE){
  G <- as.matrix(G)
  tH <- t(H)
  mtt1 <- array(0,c(T,NF))
  mtt <- array(0,c(T+1,NF))
  Rtt1 <- array(0,c(T,NF,NF))
  Rtt <- array(0,c(T+1,NF,NF))
  simAlpha <- array(0,c(T+1,NF))
  Innt <- array(0,c(T,N))
  PrecInnt <- array(0,c(T,N,N))     
  ##Forward Filtering
  Rtt[1,,] <- as.matrix(Sigma)
  for(t in 1:T){
    mtt1[t,] <- G%*%mtt[t,] 
    Rtt1[t,,] <- Sigma+G%*%Rtt[t,,]%*%t(G)
    Rtt1[t,,] <- (Rtt1[t,,]+t(Rtt1[t,,]))/2 
    TM1 <- H%*%Rtt1[t,,]
    Mti <- solve(Omega+TM1%*%tH) 
    PrecInnt[t,,] <- Mti
    TM2 <- Rtt1[t,,]%*%tH%*%Mti
    Innt[t,] <- y[t,]-lp[t,]-H%*%mtt1[t,]
    mtt[t+1,] <- mtt1[t,]+TM2%*%(Innt[t,])
    Rtt[t+1,,] <- Rtt1[t,,] - TM2%*%TM1
    Rtt[t+1,,] <- (Rtt[t+1,,]+t(Rtt[t+1,,]))/2 
  }
  if(BwSp){##Backward Simulation
    Rtt[T+1,,]=(Rtt[T+1,,]+t(Rtt[T+1,,]))/2
    simAlpha[T+1,] <- rmvnorm(1, mean = mtt[T+1,], sigma = Rtt[T+1,,],method="chol")
    for(t in T:1){
      tm <- Rtt[t,,]%*%t(G)%*%solve(Rtt1[t,,]) 
      Rt <- Rtt[t,,]-tm%*%G%*%Rtt[t,,]
      Rt <- (Rt+t(Rt))/2 
      mt <- mtt[t,]+tm%*%(simAlpha[t+1,]-mtt1[t,])
      simAlpha[t,] <- rmvnorm(1, mean = mt, sigma = Rt,method="chol")
    }
  }
  if(lglk){
    ll <- 0
    for(t in 1:T){
      ll <- ll+determinant(PrecInnt[t,,],logarithm=TRUE)$modulus[[1]]-t(Innt[t,])%*%PrecInnt[t,,]%*%Innt[t,]
    }
    ll <- ll/2-T*NF*log(2*pi)/2
  }
  ret <- list()
  if(BwSp) ret <- c(ret,list(simAlpha=simAlpha[-1,]))
  if(lglk) ret <- c(ret,list(ll=ll))
  if(filt) ret <- c(ret,list(mtt=mtt[-1,]))
  return(ret)
}



ffbs.spectral <- function(w=NULL,wFT=NULL,spec=NULL,Gvec=NULL,tau2=NULL,par=NULL,n,T,lglk=FALSE,BwSp=TRUE,NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1){
   if(is.null(wFT)){
     wFT <- real.fft.TS(TSmat.to.vect(w),n=n,T=T,inv=TRUE)
   }else{
     if(class(wFT)=="matrix"){
       wFT <- TSmat.to.vect(wFT)
     }
     if(length(wFT)!=(n*n*T)){
       print("Error: 'wFT' needs to be a vector of length T*n*n or a matrix of dimension T x n*n.")
       return()
     }
   }
   if((is.null(Gvec) & is.null(par)) | (is.null(spec) & is.null(par))){
     print("Either 'Gvec' and 'spec' and 'tau2' or, alternartively, 'par' must be supplied.")
     return()
   }
   if(is.null(Gvec)){
    wave <- wave.numbers(n)$wave
    Gvec<- get.propagator.vec(wave=wave,indCos=indCos,zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],muX=par[7],muY=par[8],dt=dt,ns=ns)
  }
   if(is.null(spec)){
     if(is.null(wave)) wave <- wave.numbers(n)$wave
     spec <- innov.spec(wave=wave,n=n,ns=ns,rho0=par[1],sigma2=par[2],zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],nu=nu,dt=dt,norm=TRUE)
   }
   if(is.null(tau2)) tau2=par[9]
   ffbsC <- .C("ffbs_spectral", wFT=as.double(wFT), bw=as.double(sum(BwSp)), ll=as.double(sum(lglk)),as.double(spec[1:ns]),as.double(Gvec$G11C),as.double(spec[indCos]),as.double(Gvec$G11),as.double(Gvec$G12),as.double(spec),as.double(tau2), as.integer(T), as.integer((NF-ns)/2), as.integer(ns))
  ret <- list()
  if(BwSp) ret <- c(ret,list(simAlpha=vect.to.TSmat(ffbsC$wFT,T=T)))
  if(lglk) ret <- c(ret,list(ll=ffbsC$ll))
  return(ret)
}


sample.four.coef <- function(w=NULL,wFT=NULL,spec=NULL,Gvec=NULL,tau2=NULL,par=NULL,n,T,NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1){
  alpha <- ffbs.spectral(w=w,wFT=wFT,spec=spec,Gvec=Gvec,tau2=tau2,par=par,n=n,T=T,lglk=FALSE,BwSp=TRUE,NF=NF,indCos=indCos,ns=ns,nu=nu,dt=dt)$simAlpha
  return(alpha)
}

loglike <- function(par=NULL,w=NULL,wFT=NULL,x=NULL,spec=NULL,Gvec=NULL,tau2=NULL,n,T,NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1,logScale=FALSE,logInd=c(1,2,3,4,5,9),negative=FALSE){
  if(logScale) par[logInd] <- exp(par[logInd])
  if(is.null(tau2)) tau2=par[9]
   if(!is.null(x)){
    lp <- apply(x,c(2,3),lin.pred,b=par[-c(1:9)])
  }else{
    lp <- 0
  }
  if(is.null(wFT)){
     wFT <- real.fft.TS(TSmat.to.vect(w-lp),n=n,T=T,inv=TRUE)
  }
   if((is.null(Gvec) & is.null(par)) | (is.null(spec) & is.null(par))){
     print("Either 'Gvec' and 'spec' and 'tau2' or, alternartively, 'par' must be supplied.")
     return()
   }
   if(is.null(Gvec)){
    wave <- wave.numbers(n)$wave
    Gvec<- get.propagator.vec(wave=wave,indCos=indCos,zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],muX=par[7],muY=par[8],dt=dt,ns=ns)
  }
   if(is.null(spec)){
     if(is.null(wave)) wave <- wave.numbers(n)$wave
     spec <- innov.spec(wave=wave,n=n,ns=ns,rho0=par[1],sigma2=par[2],zeta=par[3],rho1=par[4],gamma=par[5],alpha=par[6],nu=nu,dt=dt,norm=TRUE)
   }
  ll <- ffbs.spectral(w=w,wFT=wFT,spec=spec,Gvec=Gvec,tau2=tau2,par=par,n=n,T=T,lglk=TRUE,BwSp=FALSE,NF=NF,indCos=indCos,ns=ns,nu=nu,dt=dt)$ll
  if(negative) return(-ll) else return(ll)
}

lin.pred <- function(x,beta){
  return(x%*%beta)
}
trace.plot <- function(data,true=NULL,BurnIn=NULL,BurnInAdaptive=NULL){
  for(i in 1:dim(data)[1]){
    plot(1:dim(data)[2],data[i,],type="l",main=dimnames(data)[[1]][i])
    if(!is.null(true)){
      abline(h=true[i],lwd=2,lty=2)
    }
    if(!is.null(BurnIn)){
      abline(v=BurnIn,lwd=2,lty=2)
    }
     if(!is.null(BurnInAdaptive)){
      abline(v=BurnInAdaptive,lwd=2,lty=3)
    }
  }
}
post.dist.hist <- function(data,true=NULL,breaks=20,mean=FALSE,median=TRUE){
  for(i in 1:dim(data)[1]){
    if(!is.null(true)){
      xlim <- c(min(data[i,],true[i]),max(data[i,],true[i]))
      if(true[i]<min(data[i,])) xlim <- c(min(data[i,],true[i]-0.01*abs(true[i])),max(data[i,],true[i]))
      if(true[i]>max(data[i,])) xlim <- c(min(data[i,],true[i]),max(data[i,],true[i]+0.01*abs(true[i])))
      hist(data[i,],main=dimnames(data)[[1]][i],breaks=breaks,xlab=names(data)[i],xlim=xlim)
      abline(v=true[i],lwd=3,lty=2)
    }else{
      hist(data[i,],main=dimnames(data)[[1]][i],breaks=breaks,xlab=names(data)[i])
    }
    if(mean) abline(v=mean((data[i,])),lwd=3)
    if(median) abline(v=median((data[i,])),lwd=3)
  }
}
tobit.lambda.log.full.cond <- function(y,z,tau2,lambda){
    indNA <- is.na(y)
    ind0 <- y==0
    ind0[is.na(ind0)] <- FALSE
    indnc <- !(indNA | ind0)
  ret <- sum(log(y[indnc]^(1/lambda-1)/lambda))-sum((y[indnc]^(1/lambda)-z[indnc])^2)/2/tau2  
  return(ret)
}
mcmc.summary <- function(data,probs=c(0.025,0.5,0.975),mean=FALSE){
  ret <- matrix(0,ncol=length(probs),nrow=dim(data)[1]) 
  for(i in 1:dim(data)[1]){
    ret[i,] <- quantile(data[i,],probs=probs)
  }
  colnames(ret) <- paste(probs,rep("quant",length(probs)))
  rownames(ret) <- dimnames(data)[[1]]
  if(mean){
    means <- apply(data,1,mean)
    ret <- cbind(ret,means)
    colnames(ret) <- c(paste(probs,rep("quant",length(probs))),"mean")
  }
  return(ret)
}
print.spateMCMC <- function(x,...){
  post <- x[[1]][x$indEst,-c(1:x$BurnIn)]
  quant <- mcmc.summary(post,probs=c(0.5,0.025,0.975))
  colnames(quant) <- c("Median","2.5 %", "97.5 %")
  rownames(quant) <- dimnames(post)[[1]]
  cat("\nPosterior of parameters:\n")
  print(quant)
  cat(paste("\nResults based on ",dim(x[[1]])[2]-x$BurnIn," MCMC samples after a burn-in of ",x$BurnIn," samples\n",sep=""))
}
plot.spateMCMC <- function(x,...,trace=TRUE,hist=TRUE,medianHist=TRUE,pairs=FALSE,ask=TRUE,ToFile=FALSE,path=NULL,file=NULL,true=NULL,BurnInAdaptive=NULL,postProcess=FALSE){
  post <- x[[1]][x$indEst,]
  true <- true[x$indEst]
  if(ToFile) ask <- FALSE
  nx <- 3
  ny <- (dim(post)[1]-dim(post)[1]%%nx)/nx+1
  if(dim(post)[1]%%nx==0) ny <- ny-1
  if(trace){
    if(ToFile) jpeg(paste(path,"Traces_",file,".jpeg",sep=""),width = 750, height = 500)
    par(mfrow=c(ny,nx),oma=c(0,0,0,0),mar=c(2,2,2,0.2),ask=ask)
    trace.plot(post,BurnIn=x$BurnIn,true=true,BurnInAdaptive=BurnInAdaptive)
    if(ToFile) dev.off()
  }
  if(hist &  dim(post)[2]>(x$BurnIn+100)){
    if(ToFile) jpeg(paste(path,"PostDist_",file,".jpeg",sep=""),width = 500, height = 500)
    par(mfrow=c(ny,nx),oma=c(0,0,0,0),mar=c(2,2,2,0.2),ask=ask)
    post.dist.hist(x[[1]][,-c(1:x$BurnIn)],true=true,median=medianHist)
    if(ToFile) dev.off()
  }
  if(pairs){
    end <- dim(x[[1]])[2]
    if(end>x$BurnIn){
      spl <- sample.int(end-x$BurnIn,min(end-x$BurnIn,500))+x$BurnIn
    }else{
      spl <- sample.int(end,min(end,500))
    }
    PostSample <- x[[1]][,spl]
    col="#00009950"
    par(ask=ask)
    if(ToFile) jpeg(paste(path,"Pairs_",file,".jpeg",sep=""),width = 1000, height = 1000)
    pairs(t(PostSample),pch=".",cex=4,col=col)
    if(ToFile) dev.off()
  }
  if(postProcess){
    mean <- apply(x$xiPost,c(1,2),mean)
    spate.plot(xi=mean,ToFile=ToFile,path=path,file=file,...)
  }
}
spate.mcmc <- function(y, coord=NULL, lengthx=NULL, lengthy=NULL,Sind=NULL, n=NULL, IncidenceMat=FALSE, x = NULL, SV = c(rho0=0.2,sigma2=0.1,zeta=0.25,rho1=0.2,gamma=1,alpha=0.3,muX=0,muY=0,tau2=0.005), betaSV = rep(0, dim(x)[1]), RWCov = NULL, parh=NULL,tPred=NULL,sPred=NULL,P.rho0=Prho0,P.sigma2=Psigma2,P.zeta=Pzeta,P.rho1=Prho1,P.gamma=Pgamma,P.alpha=Palpha,P.mux=Pmux,P.muy=Pmuy,P.tau2=Ptau2, lambdaSV = 1, sdlambda = 0.01, P.lambda=Plambda, DataModel = "Normal",DimRed=FALSE,NFour=NULL, indEst = 1:9, Nmc = 10000, BurnIn = 1000, path=NULL, file=NULL, SaveToFile = FALSE, PlotToFile = FALSE, FixEffMetrop = TRUE, saveProcess = FALSE, Nsave = 200, seed = NULL, Padding = FALSE, adaptive = TRUE, NCovEst = 500, BurnInCovEst = 500, MultCov=0.5, printRWCov=FALSE, MultStdDevLambda=0.75, Separable = FALSE, Drift = !Separable, Diffusion = !Separable, logInd = c(1, 2, 3, 4, 5, 9), nu = 1, plotTrace = TRUE, plotHist = FALSE, plotPairs = FALSE, trueVal = NULL,plotObsLocations=FALSE,trace=TRUE,monitorProcess=FALSE,tProcess=NULL,sProcess=NULL){
  if(is.null(logInd)) logL="" else logL="Log_"##Update certain parameters on the log scale?
  indNN <- c(1, 2, 3, 4, 5, 9)[-match(logInd,c(1, 2, 3, 4, 5, 9))]##Which parameters cannot be negative (those not simulated on log scale)
  if(!Diffusion){
    indEst <- indEst[-match(4:6,indEst)]
    SV[4:6] <- c(0,1,0)
    logInd <- logInd[-match(4:5,logInd)]
  }
  if(!Drift){
    indEst <- indEst[-match(7:8,indEst)]
    SV[7:8] <- c(0,0)
  }
  if(is.null(RWCov)){
    noRWCov <- TRUE
    FirstFit=TRUE
  }else{
    noRWCov <- FALSE
    FirstFit=FALSE
  }
  if(noRWCov) RWCov <- diag(rep(0.005,9))
  if(is.null(Sind) & is.null(coord)){
    if((sqrt(dim(y)[2])-floor(sqrt(dim(y)[2])))!=0){
      cat("Error: the data 'y' must be on a grid. For this, the square root of dim(y)[2] must be an integer corresponding to the number of points 'n' per axis. Either put 'NA's at the points where no observations are available or, alternatively, specify the the observations locations in 'Sind' or 'coord' and specify the number of grid points 'n' per axis. \n")
      return()
    }else{
      if(is.null(n)){
        n <- sqrt(dim(y)[2])
      }else{
        if(n!=sqrt(dim(y)[2])){
          cat("Error: if the data 'y' is given on a grid, the number of grid points 'n' per axis must equal the square root of dim(y)[2]. Either put 'NA's at the points where no observations are available or, alternatively, specify the the observations locations in 'Sind' or 'coord' and specify the number of grid points 'n' per axis. \n")
          return()
        }
      }
    }
  }else{
    if(is.null(n)){
      cat("Error: the number of points per axis 'n' must be specified. \n")
      return()
    }
  }
  nOrig <- n
  T <- dim(y)[1]
  Nobs <- dim(y)[2]
  if(Padding){
    indPad <- as.vector(apply(matrix(1:n*(2*n)-2*n),1,rep,times=n))+rep(1:n,n)
    n <- 2*n
  }
  ##Padding for data on grid
  if(Padding & is.null(coord) & is.null(Sind)){
    yp <- matrix(nrow=T,ncol=n*n)
    yp[,indPad] <- y
    y <- yp
    if(!is.null(x)){
      if(dim(x)[3]!=n*n){
        cat("Error: covariates need to be available in the increased domain due to the data augmentation in the MCMC algorithm. Alternatively, use an incidence matrix approach ('IncidenceMat=TRUE', see the help) in combination with dimension reduction so that the covariates need only be available at the locations where observations are made. \n")
        return()
      }
    }
  }
  ##Observations not on grid
  if(!is.null(coord) | !is.null(Sind)){
    cellx <- rep(1:nOrig,nOrig)/nOrig-1/2/nOrig
    celly <- as.vector(apply(matrix(1:nOrig),1,rep,times=nOrig))/nOrig-1/2/nOrig 
    if(is.null(Sind)){
      if(is.null(lengthx)) coord[,1] <- coord[,1]/max(coord[,1]) else coord[,1] <- coord[,1]/lengthx
      if(is.null(lengthy)) coord[,2] <- coord[,2]/max(coord[,2]) else coord[,2] <- coord[,2]/lengthy
      Sind <- rep(0,dim(coord)[1])
      for(i in 1:dim(coord)[1]){
        d=abs(cellx-coord[i,1])^2+abs(celly-coord[i,2])^2
        Sind[i] <- which(d==min(d))
      }
    }
    if(Padding){
      SindOrig <- Sind##Indices relating the original smaller grid to the observation locations
      Sind <- indPad[Sind]##Indices relating the padded grid to the observation locations
      cellx <- rep(1:n,n)/n-1/2/n
      celly <- as.vector(apply(matrix(1:n),1,rep,times=n))/n-1/2/n
    }
    if(trace) cat("Observation locations mapped to grid. \n")
    if(plotObsLocations){
      if(PlotToFile) jpeg(paste(path,"Locations2Grid_",file,".jpeg",sep=""),width = 1000, height = 1000) else par(ask=T)
      plot(cellx,celly,pch=".",cex=1,xlab="x",ylab="y",xlim=c(0,1),ylim=c(0,1),main="Grid cell and observation points")
      abline(v=(0:n)/n)
      abline(h=(0:n)/n)
      points(cellx[Sind],celly[Sind],pch=".",cex=5)
      if(!is.null(coord)) points(coord[,1],coord[,2],pch=18)
      if(PlotToFile) dev.off()
    }
  }
  if(!is.null(NFour)){
    if(NFour<n*n & !DimRed){
      DimRed <- TRUE
      cat("Warning: 'DimRed' was changed to 'TRUE' since NFour<n*n implies a dimension reduction. \n") 
    }
  }
  if(DimRed) NF <- NFour else NF <- n*n
  ##SPDE initialization
  spateFT <- spate.init(n=n,T=T,NF=NF)
  if(IncidenceMat){
    if(NFour==n*n | is.null(NFour)){
      cat("Error: If 'IncidenceMat' equals TRUE dimension reduction needs to be done since the FFT cannot be used. Please specify 'NFour' and make sure that NFour<<n*n. \n")
      return()
    }
    Phi <- get.real.dft.mat(wave=spateFT$wave, indCos=spateFT$indCos, ns=spateFT$ns, n=n)
    I <- matrix(0,length(Sind),dim(Phi)[1])
    for(i in 1:length(Sind)){
      I[i,Sind[i]] <- 1
    }
    IPhi <- I%*%Phi
    rm(I)
    if(trace) cat("Incidence matrix constructed. \n")
  }
  ##Data
  wh <- y
  indNA <- is.na(y)
  ParNames <- c("rho_0","sigma^2","zeta","rho_1","gamma","alpha","mu_x","mu_y","tau^2")
  if(is.null(parh)){
    parh <- array(0,c(length(ParNames),Nmc+1))
    parh[,1] <- SV
    dimnames(parh)[[1]] <- ParNames
    if(!FixEffMetrop) dimnames(RWCov)[[2]] <- dimnames(RWCov)[[1]] <- ParNames
    Prediction=FALSE
  }else{
    Prediction=TRUE
    if(!is.null(x)) indFECoef <- (length(ParNames)+1):(length(ParNames)+dim(x)[1])
    ypred <- array(0,c(length(tPred),length(sPred),Nmc-BurnIn))##change indS with padding
  }
  ##Covariates
  if(!is.null(x) & !Prediction){
    Ncov <- dim(x)[1]
    CovNames <- dimnames(x)[[1]]
    if(is.null(CovNames)) CovNames <- paste(rep("Cov",Ncov),1:Ncov)
    betah <- matrix(0,nrow=Ncov,ncol=Nmc+1)
    betah[,1] <- betaSV
    lp <- apply(x,c(2,3),lin.pred,b=betaSV)
    if(FixEffMetrop){##Covariates if the coefficients of the fixed effects are updated jointly in te MH step
      parh <- array(0,c(length(ParNames)+dim(x)[1],Nmc+1))
      parh[,1] <- c(SV,betaSV)
      dimnames(parh)[[1]] <- c(ParNames,CovNames)
      indFECoef <- (length(ParNames)+1):(length(ParNames)+dim(x)[1])
      indEst <- c(indEst,indFECoef)
      if(noRWCov) RWCov=diag(c(diag(RWCov),rep(0.001,Ncov)))
      dimnames(RWCov)[[2]] <- dimnames(RWCov)[[1]] <- c(ParNames,CovNames)
    }else{##Covariates if Gibbs sampling is used
      betah[,1] <- betaSV
      xv <- apply(x,1,TSmat.to.vect)
      xTxi <- solve(crossprod(xv))
    }
  }else{##No covariates
    lp <- lpV <- array(0,c(dim(wh)[1],dim(wh)[2]))
  }
  ##Processes
  if(IncidenceMat) xih <- array(0,c(T,Nobs)) else xih <- array(0,c(T,n*n))
  alphah <- array(0,c(T,NF))
  xiPost <- NULL
  if(saveProcess){
    if(Nsave*T*n*n*8/1000000>100) cat(paste("Warning: saving ",Nsave," samples from the posterior of the process requires ",Nsave*T*n*n*8/1000000," MB of memory."," \n",sep=""))
    xiPost <-array(0,c(T,n*n,Nsave))
  }
  if(DataModel=="SkewTobit"){
    ind0 <- y==0
    ind0[is.na(ind0)] <- FALSE
    indnc <- !(indNA | ind0)
    if(!Prediction){
      wh[indnc] <- y[indnc]^(1/lambdaSV)
      lambdah <- rep(lambdaSV,Nmc+1)
    }
  }
  ##Acceptance rate
  AcRate <-  AcRate2 <- rep(0,1)
  if(DataModel=="SkewTobit") AcRate <-  AcRate2 <- rep(0,2)
  if(!is.null(seed)) set.seed(seed)
  iFirst <- TRUE
  if(monitorProcess){
    if(is.null(tProcess) | is.null(sProcess)){
      cat("Error: if 'monitorProcess=TRUE' is selected, 'tProcess' and 'sProcess' (temporal and spatial locations) for monitoring xi need to be specified. \n")
      return()
    }
    Ntrace <- max(length(tProcess),length(sProcess))
    if(length(tProcess)<Ntrace){
      tProcess <- c(tProcess,sample.int(dim(xih)[1],Ntrace-length(tProcess)))
      cat("Warning: The indeces 'tProcess' and 'sProcess' (temporal and spatial locations) for monitoring xi in the MCMC algorithm do not have the same length. Random numbers have been added to the shorter index so that they have the same length. \n")
    }
    if(length(sProcess)<Ntrace){
      sProcess <- c(sProcess,sample.int(dim(xih)[2],Ntrace-length(sProcess)))
      cat("Warning: The indeces 'tProcess' and 'sProcess' (temporal and spatial locations) for monitoring xi in the MCMC algorithm do not have the same length. Random numbers have been added to the shorter index so that they have the same length. \n")
    }
    xiTrace  <-array(0,c(Ntrace,Nmc+1))
    indXiTrace <- array(0,c(2,Ntrace))
    indXiTrace[1,] <- tProcess
    indXiTrace[2,] <- sProcess
    dimnames(xiTrace)[[1]] <- paste("Trace plot of xi at time ",tProcess," and location ",sProcess,sep="")
  }
  if(trace) cat("Starting MCMC algorithm: \n")
  
  ##Start MCMC algorithm
  for(i in 1:Nmc){
    if(Prediction){
      if(DataModel=="SkewTobit") wh[indnc] <- y[indnc]^(1/parh["lambda",i])
      if(!is.null(x)){
        if(IncidenceMat & Nobs!=dim(x)[2]) lp <- apply(x[,,SindOrig],c(2,3),lin.pred,b=parh[indFECoef,i]) else lp <- apply(x,c(2,3),lin.pred,b=parh[indFECoef,i])
      }
    }
    ##Latent variable Z
    zh <- lp+xih 
    ##Gibbs step for NAs
    if(sum(indNA)>0) wh[indNA] <- rnorm(sum(indNA),mean=zh[indNA],sd=sqrt(parh[9,i]))
    ##Tobit model: M-H step for transformation parameter
    if(DataModel=="SkewTobit" & !Prediction){
      wh[ind0] <- rtruncnorm(sum(ind0),mean=zh[ind0],sd=sqrt(parh[9,i]),b=rep(0,sum(ind0)))
      lambdaV <- exp(log(lambdah[i])+rnorm(1,mean=0,sd=sdlambda))
      al <- min(1, exp(tobit.lambda.log.full.cond(y=y,z=zh,tau2=parh[9,i],lambda=lambdaV)-tobit.lambda.log.full.cond(y=y,z=zh,tau2=parh[9,i],lambda=lambdah[i])+P.lambda(lambdaV,log=TRUE)-P.lambda(lambdah[i],log=TRUE))*lambdaV/lambdah[i])
      if(runif(1)<=al){
        lambdah[i+1] <- lambdaV
        wh[indnc] <- y[indnc]^(1/lambdah[i+1])
        AcRate[2] <- AcRate[2]+1
        if(i>(BurnInCovEst+NCovEst)) AcRate2[2] <- AcRate2[2]+1
      }else{
        lambdah[i+1] <- lambdah[i]
      }
    }
    ##Gibbs step for coefficients beta
    if(!is.null(x) & !FixEffMetrop & !Prediction){
      mbeta <- xTxi%*%t(xv)%*%TSmat.to.vect(wh-xih)
      betah[,i+1] <- rmvnorm(1, mean = mbeta, sigma = parh[9,i]*xTxi)
      lp <- apply(x,c(2,3),lin.pred,b=betah[,i+1])
    }
    ##Initialization of 'xih' or making predictions when in prediction mode
    if(i==1 | Prediction){
      spec <- innov.spec(wave=spateFT$wave,n=n,ns=spateFT$ns,rho0=parh[1,i],sigma2=parh[2,i],zeta=parh[3,i],rho1=parh[4,i],gamma=parh[5,i],alpha=parh[6,i],nu=nu,dt=1,norm=TRUE)
      if(!IncidenceMat){
        if(DimRed) wFT <- TSmat.to.vect(vect.to.TSmat(real.fft.TS(TSmat.to.vect(wh-lp),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT),T=T)[spateFT$IndFour,]) else wFT <- real.fft.TS(TSmat.to.vect(wh-lp),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT)
        Gvec <- get.propagator.vec(wave=spateFT$wave,indCos=spateFT$indCos,zeta=parh[3,i],rho1=parh[4,i],gamma=parh[5,i],alpha=parh[6,i],muX=parh[7,i],muY=parh[8,i],dt=1,ns=spateFT$ns)
        alphah <- ffbs.spectral(wFT=wFT,spec=spec,Gvec=Gvec,tau2=parh[9,i],n=n,T=T,lglk=FALSE,BwSp=TRUE,NF=dim(spateFT$wave)[2],indCos=spateFT$indCos,ns=spateFT$ns)$simAlpha
        if(DimRed){
          alphahC <- array(0,c(T,n*n))
          alphahC[,spateFT$IndFour] <- alphah
          xih <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphahC),n=n,T=T,inv=FALSE,indFFT=spateFT$indFFT),T=T)
        }else{
          xih <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphah),n=n,T=T,inv=FALSE,indFFT=spateFT$indFFT),T=T)
        }
      }else{##IncidenceMat
        G <- get.propagator(wave=spateFT$wave,indCos=spateFT$indCos,zeta=parh[3,i],rho1=parh[4,i],gamma=parh[5,i],alpha=parh[6,i],muX=parh[7,i],muY=parh[8,i],dt=1,ns=spateFT$ns)
        alphah <- ffbs(y=wh,lp=lp,G=G,Sigma=diag(spec),H=IPhi,Omega=diag(rep(parh[9],dim(wh)[2])), lglk=FALSE,BwSp=TRUE)$simAlpha
        xih=t(IPhi%*%t(alphah))
      }
      if(Prediction){
        if(IncidenceMat & length(sPred)!=Nobs){
          alphahC <- array(0,c(T,n*n))
          alphahC[,spateFT$IndFour] <- alphah
          xihPred <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphahC),n=n,T=T,inv=FALSE,indFFT=spateFT$indFFT),T=T)
          xihPred <- xihPred[,indPad[sPred]]
          lpPred <- apply(x,c(2,3),lin.pred,b=parh[indFECoef,i])
        }else{
          xihPred <- xih
          lpPred <- lp
        }
        zh <- lpPred+xihPred 
        if(i>BurnIn){
          ypred[,,i-BurnIn] <- zh[tPred,sPred]
          if(DataModel=="SkewTobit"){
            ypred[,,i-BurnIn][ypred[,,i-BurnIn]<0] <- 0
            ypred[,,i-BurnIn] <- ypred[,,i-BurnIn]^parh["lambda",i]
          }
        }
      }
    }
    
    if(!Prediction){##M-H step for joint Update (JU) of process 'xih' and parameters
      m=parh[,i]
      m[logInd]=log(parh[logInd,i])
      parV <- m
      parV[indEst] <- rmvnorm(1, mean = m[indEst], sigma = as.matrix(RWCov[indEst,indEst]), method="chol")
      if(logL=="Log_") parV[logInd]=exp(parV[logInd])
      if(sum(parV[indNN]<=0)>0| parV[5]<1e-6 | parV[5]>1e6 | parV[3]<1e-8){
        al <- 0
      }else{
        specV <- innov.spec(wave=spateFT$wave,n=n,ns=spateFT$ns,rho0=parV[1],sigma2=parV[2],zeta=parV[3],rho1=parV[4],gamma=parV[5],alpha=parV[6],nu=nu,dt=1,norm=TRUE)
        if(!IncidenceMat){
          if(DimRed) wFT <- wFTV <- TSmat.to.vect(vect.to.TSmat(real.fft.TS(TSmat.to.vect(wh-lp),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT),T=T)[spateFT$IndFour,]) else wFT <- wFTV <- real.fft.TS(TSmat.to.vect(wh-lp),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT)
          if(FixEffMetrop & !is.null(x)){
            lpV <- apply(x,c(2,3),lin.pred,b=parV[indFECoef])
            if(DimRed) wFTV <- TSmat.to.vect(vect.to.TSmat(real.fft.TS(TSmat.to.vect(wh-lp),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT),T=T)[spateFT$IndFour,]) else wFTV <- real.fft.TS(TSmat.to.vect(wh-lpV),n=n,T=T,inv=TRUE,indFFT=spateFT$indFFT)
          }
          GvecV <- get.propagator.vec(wave=spateFT$wave,indCos=spateFT$indCos,zeta=parV[3],rho1=parV[4],gamma=parV[5],alpha=parV[6],muX=parV[7],muY=parV[8],dt=1,ns=spateFT$ns)
          mllV <- ffbs.spectral(wFT=wFTV,spec=specV,Gvec=GvecV,tau2=parV[9],n=n,T=T,lglk=TRUE,BwSp=FALSE,NF=dim(spateFT$wave)[2],indCos=spateFT$indCos,ns=spateFT$ns)$ll
          mll <- ffbs.spectral(wFT=wFT,spec=spec,Gvec=Gvec,tau2=parh[9,i],n=n,T=T, lglk=TRUE,BwSp=FALSE,NF=dim(spateFT$wave)[2],indCos=spateFT$indCos,ns=spateFT$ns)$ll
        }else{##IncidenceMat
          GV <- get.propagator(wave=spateFT$wave,indCos=spateFT$indCos,zeta=parV[3],rho1=parV[4],gamma=parV[5],alpha=parV[6],muX=parV[7],muY=parV[8],dt=1,ns=spateFT$ns)
          if(FixEffMetrop & !is.null(x)) lpV <- apply(x,c(2,3),lin.pred,b=parV[indFECoef]) else lpV <- lp
          ffbs <- ffbs(y=wh,lp=lpV,G=GV,Sigma=diag(specV),H=IPhi,Omega=diag(rep(parV[9],dim(wh)[2])), lglk=TRUE,BwSp=TRUE)
          mllV <- ffbs$ll
          mll <- ffbs(y=wh,lp=lp,G=G,Sigma=diag(spec),H=IPhi,Omega=diag(rep(parh[9,i],dim(wh)[2])), lglk=TRUE,BwSp=FALSE)$ll
        }
        if(sum(4==indEst)>0) al <- min(1, exp(P.rho0(parV[1],log=TRUE)+P.sigma2(parV[2],log=TRUE)+P.zeta(parV[3],log=TRUE)+P.rho1(parV[4],log=TRUE)+P.gamma(parV[5],log=TRUE)+P.alpha(parV[6],log=TRUE)+P.mux(parV[7],log=TRUE)+P.muy(parV[8],log=TRUE)+P.tau2(parV[9],log=TRUE)-(P.rho0(parh[1,i],log=TRUE)+P.sigma2(parh[2,i],log=TRUE)+P.zeta(parh[3,i],log=TRUE)+P.rho1(parh[4,i],log=TRUE)+P.gamma(parh[5,i],log=TRUE)+P.alpha(parh[6,i],log=TRUE)+P.mux(parh[7,i],log=TRUE)+P.muy(parh[8,i],log=TRUE)+P.tau2(parh[9,i],log=TRUE))+mllV-mll+sum(log(parV[logInd]))-sum(log(parh[logInd,i])))) else al <- min(1, exp(P.rho0(parV[1],log=TRUE)+P.sigma2(parV[2],log=TRUE)+P.zeta(parV[3],log=TRUE)+P.gamma(parV[5],log=TRUE)+P.alpha(parV[6],log=TRUE)+P.mux(parV[7],log=TRUE)+P.muy(parV[8],log=TRUE)+P.tau2(parV[9],log=TRUE)-(P.rho0(parh[1,i],log=TRUE)+P.sigma2(parh[2,i],log=TRUE)+P.zeta(parh[3,i],log=TRUE)+P.gamma(parh[5,i],log=TRUE)+P.alpha(parh[6,i],log=TRUE)+P.mux(parh[7,i],log=TRUE)+P.muy(parh[8,i],log=TRUE)+P.tau2(parh[9,i],log=TRUE))+mllV-mll+sum(log(parV[logInd[-4]]))-sum(log(parh[logInd[-4],i]))))
      }
      if(runif(1)<=al){
        parh[,i+1] <- parV
        spec <- specV
        if(!IncidenceMat){
          Gvec <- GvecV
          alphah <- ffbs.spectral(wFT=wFTV,spec=specV,Gvec=GvecV,tau2=parV[9],n=n,T=T,lglk=FALSE,BwSp=TRUE,NF=dim(spateFT$wave)[2],indCos=spateFT$indCos,ns=spateFT$ns)$simAlpha
          if(DimRed){
            alphahC <- array(0,c(T,n*n))
            alphahC[,spateFT$IndFour] <- alphah
            xih <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphahC),n=n,T=T,inv=FALSE,indFFT=spateFT$indFFT),T=T)
          }else{
            xih <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphah),n=n,T=T,inv=FALSE,indFFT=spateFT$indFFT),T=T)
          }
        }else{##IncidenceMat
          G <- GV
          alphah <- ffbs$simAlpha
          xih=t(IPhi%*%t(alphah))
        }
        if(FixEffMetrop) lp <- lpV
        AcRate[1] <- AcRate[1]+1
        if(i>(BurnInCovEst+NCovEst)){
          AcRate2[1] <- AcRate2[1]+1
        }
      }else{
        parh[,i+1] <- parh[,i]
      }
      if(monitorProcess) xiTrace[,i+1] <- xih[t(indXiTrace)]
      if(saveProcess & i>BurnIn & (i-BurnIn)%%trunc((Nmc-BurnIn)/Nsave)==0 & (i-BurnIn)/trunc((Nmc-BurnIn)/Nsave)<=Nsave){
        if(IncidenceMat) xiPost[,,(i-BurnIn)/trunc((Nmc-BurnIn)/Nsave)] <- Phi%*%t(alphah) else xiPost[,,(i-BurnIn)/trunc((Nmc-BurnIn)/Nsave)] <- xih
      }
      ##Proposal covariance estimation
      if((i%%100==0) & AcRate[1]<=0.01){
        cat(paste("Acceptance rate for hyperparameters random walk is less than 1% after",i,"iterations. Because of this, the proposal covariance matrix is devided by 100. \n"))
        RWCov <- RWCov/100
      }
      if(i==200 & AcRate[1]==0) RWCov <- RWCov/100
      if(i>=(BurnInCovEst+NCovEst) & (i-BurnInCovEst)%%NCovEst==0 & adaptive){
        if(logL=="Log_"){
          parhC <- parh
          parhC[logInd,] <- log(parhC[logInd,])
        }
        if(i>=(BurnInCovEst+NCovEst)){
          RWCovP <- MultCov*cov(t(parhC[,(BurnInCovEst+1):(i+1)]))
          eigenv <- eigen(RWCovP[indEst,indEst])$value
          if(is.double(eigenv) &  sum(eigenv<=0)==0){
            RWCov <- RWCovP
            if(printRWCov){
              cat("Estimated proposal covariance for hyperparameters: \n")
              print(signif(RWCov[indEst,indEst],digits=2))
            }
          }else{
            cat("Warning: random walk step for hyperparameters: estimated proposal covariance matrix is not positive definite. MCMC algorithm is continued with current proposal covariance matrix. Alternatively restart the algorithm with different settings (initial values, initial proposal covariance matrix etc.). \n")
          }
        }
        if(DataModel=="SkewTobit"){
          sdamma <- MultStdDevLambda*sd(lambdah[(BurnInCovEst+1):(i+1)])
        }
      }
    }
    if(i==1) t1 <- Sys.time()
    if(trace & ((i%%10==0 & i<=400) | (i%%20==0 & i>400 & i<=2000)| (i%%50==0 & i>2000 & i<=10000)| (i%%100==0 & i>10000 & i<=20000)| (i%%500==0 & i>20000)| i==Nmc)) cat(".")
    if(i==50 | (i%%100==0 & i<=400) | (i%%200==0 & i>400 & i<=2000)|(i%%500==0 & i>2000 & i<=10000)| (i%%1000==0 & i>10000 & i<=20000)| (i%%5000==0 & i>20000)| i==Nmc){
      if(trace){
        cat(paste("Iteration number: ",i," \n",sep=""))
        t2 <- Sys.time()
        dt <- (t2-t1)/i*(Nmc-i)
        units(dt) <- "secs"
        if(dt<=60){
          tdif <- paste(round(dt,digits=3),"secs")
        }else if(dt>60 & dt<=(60^2)){
          units(dt) <- "mins"
          tdif <- paste(round(dt,digits=3),"mins")
        }else if(dt>(60^2) & dt<=(60^2*24)){
          units(dt) <- "hours"
          tdif <- paste(round(dt,digits=3),"hours")
        }else{
          units(dt) <- "days"
          tdif <- paste(round(dt,digits=3),"days")
        }
        dti <- (t2-t1)/i
        units(dti) <- "secs"
        cat(paste("COMPUTING TIME for one iteration: ",round(dti,digits=3)," secs. Estimated remaining computing time: ",tdif," \n",sep=""))
        if(!Prediction){
          if(i<=(BurnInCovEst+NCovEst)){
            cat("ACCEPTANCE RATE for Metropolis-Hastings step:")
            cat(paste("  Hyperparameters: ",round(AcRate[1]/i,digits=2),sep=""))
            if(DataModel=="SkewTobit") cat(paste(",  Tobit transformation parameter: ",round(AcRate[2]/i,digits=2),sep=""))
            cat("\n")
          }else{
            cat("ACCEPTANCE RATE for Metropolis-Hastings step after burn-in:")
            cat(paste("  Hyperparameters: ",round(AcRate2[1]/i,digits=2),sep=""))
            if(DataModel=="SkewTobit") cat(paste(",  Tobit transformation parameter: ",round(AcRate2[2]/i,digits=2),sep=""))
            cat("\n")
          }
        }
      }
      if(!Prediction){
        SimDataRes <- parh[1:length(ParNames),1:(i+1)]
        rownames(SimDataRes) <- ParNames
        if(!is.null(x)){
          if(FixEffMetrop) betah[,1:(i+1)] <- parh[indFECoef,1:(i+1)]
          names <- rownames(SimDataRes)
          SimDataRes <- rbind(SimDataRes,betah[,1:(i+1)])
          rownames(SimDataRes) <- c(names,CovNames)
        }
        if(DataModel=="SkewTobit"){
          names <- rownames(SimDataRes)
          SimDataRes <- rbind(SimDataRes,lambdah[1:(i+1)])
          rownames(SimDataRes) <- c(names,"lambda")
        }
        indEst2=indEst
        if(!is.null(x) & !FixEffMetrop) indEst2 <- c(indEst2,10:(9+length(CovNames)))
        if(DataModel=="SkewTobit") indEst2 <- c(indEst2,dim(SimDataRes)[1])
        spateMCMC <- list(Post=SimDataRes[,-1],xiPost=xiPost,RWCov=RWCov,BurnIn=BurnIn,Padding=Padding,indEst=indEst2,nu=nu,DataModel=DataModel)
        class(spateMCMC) <- "spateMCMC"
        if(monitorProcess & !PlotToFile & !iFirst) dev.set(devPar)
        plot(spateMCMC,ToFile=PlotToFile,path=path,file=file,pairs=plotPairs,hist=plotHist,trace=plotTrace,ask=FALSE,true=trueVal,BurnInAdaptive=(BurnInCovEst+NCovEst))
        if(SaveToFile){
          save(spateMCMC,file=paste(path,file,sep="")) ##cat("spateMCMC object saved \n")
        }
        if(monitorProcess){
          nx <- ceiling(sqrt(length(sProcess)+length(tProcess)+Ntrace))
          ny <- floor(sqrt(length(sProcess)+length(tProcess)+Ntrace))
          if(PlotToFile){
            jpeg(paste(path, "ProcessSample_", file, ".jpeg", sep = ""), width = nx/ny*500, height = 500)
          }else{
            if(iFirst){
              devPar <- dev.cur()
              dev.new()
              devProc <- dev.cur()
            }else{
              dev.set(devProc)
            }
          }
          par(mfrow = c(ny, nx), oma = c(0, 0, 0, 0), mar = c(2, 2, 2, 0.2))
          trace.plot(xiTrace[,1:(i+1)])
          if(!is.null(sProcess)){
            for(s in sProcess){
              if(!IncidenceMat) main=paste("One sample of w and xi at grid point ",s,sep="") else main=paste("One sample of w and xi at station ",s,sep="")
              plot(1:min(T,150),wh[1:min(T,150),s],type="l",xlab="Time",main=main)
              abline(h=0,lty=3)
              lines(xih[1:min(T,150),s],lty=2)
              abline(v=tProcess,lwd=1,col="red")
            }
          }
          if(!is.null(tProcess)){
            for(t in tProcess){
              if(!IncidenceMat) image(1:n,1:n,matrix(xih[t,],nrow=n),main=paste("One sample of xi at time ",t,sep=""),xaxt='n',yaxt='n',col = cols()) else image(1:n,1:n,matrix(t(Phi%*%t(alphah))[t,],nrow=n),main=paste("One sample of xi at time ",t,sep=""),xaxt='n',yaxt='n',col = cols()) 
            }
          }
          if(PlotToFile) dev.off()
        }
      }
      if(iFirst) iFirst <- FALSE
    }
  }
  if(Prediction) return(ypred) else return(spateMCMC)
}


spate.predict <- function(y, tPred,sPred=NULL,xPred=NULL,yPred=NULL,spateMCMC,Nsim=200,BurnIn=5, coord=NULL, lengthx=NULL, lengthy=NULL,Sind=NULL, n=NULL, IncidenceMat=FALSE, x = NULL, DataModel = "Normal",DimRed=FALSE,NFour=NULL, seed = NULL, nu = 1,trace=FALSE){
  if(max(tPred)>dim(y)[1]) y <- rbind(y,matrix(ncol=dim(y)[2],nrow=(max(tPred)-dim(y)[1])))
  if(is.null(sPred)){
    if(is.null(xPred)){
      if(is.null(n)) n <- sqrt(dim(y)[2])
      sPred <- 1:(n^2)
    }else{
      sPred <- rep(0,length(xPred))
      for(i in 1:length(xPred)) sPred[i] <- xPred[i]+n*(yPred[i]-1)
    }
  }
  spl <- ((1:(Nsim+BurnIn))-1)*trunc((dim(spateMCMC$Post)[2]-spateMCMC$BurnIn)/(Nsim+BurnIn))+1+spateMCMC$BurnIn
  post <- spateMCMC[[1]][,spl]
  ypred <- spate.mcmc(y=y, tPred=tPred,sPred=sPred,Nmc=(Nsim+BurnIn),BurnIn=BurnIn, coord=coord, lengthx=lengthx, lengthy=lengthy,Sind=Sind, n=n, IncidenceMat=IncidenceMat, x = x, DataModel = DataModel,DimRed=DimRed,NFour=NFour, seed = seed, Padding = spateMCMC$Padding, nu = nu,trace=trace,parh=post)
  if(is.null(xPred)) return(ypred) else return(apply(ypred,3,diag))
}

