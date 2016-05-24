#Maximum likelihood estimation
#for
#t,x,n,Zb,Zw,data,erho,esigma,ebeta,ealphab,ealphaw,truth,precision
#see documentation

#@Rfun specifies function used to calculate R in the likelihood

  ei <- function(formula, total=NULL, Zb=1,Zw=1, id=NA, data=NA, erho=.5, esigma=.5,
               ebeta=.5, ealphab=NA, ealphaw=NA, truth=NA, simulate=TRUE,covariate=NULL, lambda1=4, lambda2=2, covariate.prior.list=NULL, tune.list=NULL, start.list=NULL, sample=1000, thin=1, burnin=1000, verbose=0, ret.beta="r", ret.mcmc=TRUE, usrfun=NULL){
    #Extract formula
    dv <- terms.formula(formula)[[2]]
    iv <- terms.formula(formula)[[3]]
    t <- as.character(dv)
    x <- as.character(iv)
    n <- as.character(total)
    id <- as.character(id)
    
    if(length(dv)==1){
    print("Running 2x2 ei")

    if(simulate==FALSE){
    dbuf <- ei.estimate(t,x, n,id=id,data=data, Zb=Zb, Zw=Zw, erho=erho, esigma=esigma, ebeta=ebeta, ealphab=ealphab, ealphaw=ealphaw, truth=truth)
    return(dbuf)
}
    if(simulate==TRUE){
    #If the table is two by two, use ei
    dbuf <- tryCatch(tryCatch(ei.estimate(t,x, n,id=id,
                        data=data, Zb=Zb, Zw=Zw, erho=erho,
                                          esigma=esigma, ebeta=ebeta,
                                          ealphab=ealphab,
                                          ealphaw=ealphaw,
                                          truth=truth),
                              error=function(x) ei(t,x,n, id=id,
                                      data=data, Zb=Zb, Zw=Zw,erho=3,esigma=esigma, ebeta=ebeta, ealphab=ealphab, ealphaw=ealphaw, truth=truth)), error=function(x) ei.estimate(t,x,n,
                                     id=id, data=data, Zb=Zb, Zw=Zw, erho=5,esigma=esigma, ebeta=ebeta, ealphab=ealphab, ealphaw=ealphaw, truth=truth))
    dbuf.sim <- ei.sim(dbuf)
    return(dbuf.sim)
}
  }

    if(length(dv)>1){
    print("Running eiRxC")
    #If the table is RxC use eiRxC
    dbuf <- ei.MD.bayes(formula, data=data, total=total, covariate=covariate, lambda1=lambda1, lambda2=lambda2, covariate.prior.list=covariate.prior.list, tune.list=tune.list, start.list=start.list, sample=sample, thin=thin, burnin=burnin, verbose=verbose, ret.beta=ret.beta, ret.mcmc=ret.mcmc, usrfun=usrfun)
    dbuf$data <- data
    dbuf$total <- n
    dbuf$formula <- formula
    class(dbuf) <- "ei"
    return(dbuf)
  }

}

ei.estimate <- function(t,x,n,id,Zb=1,Zw=1, data=NA, erho=.5, esigma=.5, ebeta=.5,
               ealphab=NA, ealphaw=NA, truth=NA, Rfun=2, precision=4){

#Check to make sure data is not null
  if(!missing(data)){
    t <- data[[t]]
    x <- data[[x]]
    n <- data[[n]]
    if(is.character(Zb)) Zb <- data[[Zb]]
    if(is.character(Zw)) Zw <- data[[Zw]]
    id <- data[[id]]
  }
  
  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)

  #If there are no covariates, run simplest R function
  if(dim(Zb)[1] == 1 & Zb[1,1] == 1 & dim(Zw)[1] == 1 & Zw[1,1] == 1) Rfun=5
  if (dim(Zb)[1]==1 & Zb[1,1]==1) Zb <- as.matrix(rep(1,length(x)))
  if (dim(Zw)[1]==1 & Zw[1,1]==1) Zw <- as.matrix(rep(1,length(x)))

  #Extract the number of covariates
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]

  #Starting values
  start <- c(0,0,-1.2,-1.2, 0, rep(0, numb+numw))

  message("Maximizing likelihood")
 solution <- ucminf(start, like, y=t, x=x, n=n, Zb=Zb,
     Zw=Zw,numb=numb, erho=erho, esigma=esigma,
     ebeta=ebeta, ealphab =ealphab,ealphaw=ealphaw, Rfun=Rfun, hessian=3)
#This didn't work
  #solution <- optim(start, like, y=t, x=x, n=n, Zb=Zb,
  #Zw=Zw,numb=numb, erho=erho, esigma=esigma,
  #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, hessian=T,
  #Rfun=Rfun, method="BFGS")
  #control=list(maxeval=10))
  #print(solution$par)
  #print(solution$convergence)
  #solution <- genoud(like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb,
  #erho=erho, esigma=esigma,
  #ebeta=ebeta, ealphab
  #=ealphab, ealphaw=ealphaw, nvars=5, starting.values=start)
  #solution <- maxLik(like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb,
  #erho=erho, esigma=esigma,
  #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw,start=start)
  #solution <- subplex(start, like, y=t, x=x, n=n, Zb=Zb,
  #Zw=Zw,numb=numb, erho=erho,esigma=esigma,
  #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw)
  #solution <- nlminb(start, like,y=t, x=x, n=n, Zb=Zb,
  #Zw=Zw,numb=numb, erho=erho, esigma=esigma,
  #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, Rfun=Rfun,
  #hessian=T)
  #

  #Find values of the Hessian that are 0 or 1.
  covs <- as.logical(ifelse(diag(solution$hessian)==0|
                            diag(solution$hessian)==1,0,1))
  hessian <- solution$hessian[covs,covs]
  output <- list(solution$par, solution$hessian,hessian, erho, esigma,
                 ebeta, ealphab, ealphaw, numb, x, t, n, Zb, Zw,
                 truth,precision, covs, Rfun, id)

  names(output) <- c("phi", "hessian","hessianC",  "erho",
                     "esigma", "ebeta", "ealphab", "ealphaw", "numb",
                     "x", "t", "n", "Zb", "Zw",
                     "truth", "precision", "covs", "Rfun", "id")
  class(output) <- "ei"
  return(output)
}

 ei.sim <- function(ei.object){
   hessian <- ei.object$hessianC
   erho <- ei.object$erho
   esigma <- ei.object$esigma
   ebeta <- ei.object$ebeta
   ealphab <- ei.object$ealphab
   ealphaw <- ei.object$ealphaw
   numb <- ei.object$numb
   covs <- ei.object$covs
   Rfun <- ei.object$Rfun
   x <- ei.object$x
   t <- ei.object$t
   n <- ei.object$n
   Zb <- ei.object$Zb
   Zw <- ei.object$Zw
   truth <- ei.object$truth
   id <- ei.object$id
   precision <- ei.object$precision
  #Begin Importance Sampling
  message("Importance Sampling..")
  keep <- matrix(data=NA, ncol=(length(ei.object$phi)))
  resamp <- 0
  while(dim(keep)[1] < 100){
    keep <- .samp(t,x,n, Zb, Zw, ei.object$phi, hessian, 100, keep,
                 numb=numb, covs, erho, esigma,
                 ebeta, ealphab, ealphaw, Rfun)
    resamp = resamp + 1
  }

  #Extract values from importance sampling
  keep <- keep[2:100,]
  mu <- keep[,1:2]
  sd <- keep[,3:4]
  rho <- keep[,5]
  Bb0v <- keep[,6:(5+numb)]
  Bw0v <- keep[,(6+numb):length(ei.object$phi)]
  sd[,1] <- exp(sd[,1])
  sd[,2] <- exp(sd[,2])

  #Reparamterize
  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)
  Bb0v <- as.matrix(Bb0v)
  Bw0v <- as.matrix(Bw0v)
  mu1 <- mu[,1]*(.25 + sd[,1]^2)+ .5 + t(as.matrix(apply(Zb,2,
    function (x) x - mean(x)))%*%t(Bb0v))
  mu2 <- mu[,2]*(.25 + sd[,2]^2) + .5 + t(as.matrix(apply(Zw,2,
    function (x) x - mean(x)))%*%t(Bw0v))

  #phin <- dmvnorm(psi, par, log=T)
  rho <- (exp(2*rho)-1)/(exp(2*rho) +1)
  psi <- cbind(mu1, mu2, sd, rho)
  bb <- psi[,1:length(x)]
  bw <- psi[,(length(x)+1):(length(x)*2)]
  sb <- psi[,(length(x)*2+1)]
  sw <- psi[,(length(x)*2+2)]
  rho <- psi[,(length(x)*2+3)]
  omx <- 1-x
  sbw <- rho*sb*sw
  betab <- matrix(nrow=length(x),ncol=dim(keep)[1])
  betaw <- matrix(nrow=length(x),ncol=dim(keep)[1])
  homoindx <- ifelse(x==0, 1, 0)
  homoindx <- ifelse(x==1, 2, homoindx)
  enumtol=.0001
  cT0 <- t < enumtol & homoindx==0
  cT1 <- t > (1-enumtol) & homoindx==0
  ok <- ifelse(homoindx == 0 & cT0 == 0 & cT1 == 0,T, F)
  wh <- homoindx==1
  bl <- homoindx==2
  for (i in 1:dim(keep)[1]){
    sig2 <- sb[i]^2*x^2 + sw[i]^2*omx^2 + sbw[i]*2*x*omx
    omega <- sb[i]^2*x + sbw[i]*omx
    eps <- t - (bb[i,])*x - (bw[i,])*omx
    mbb <- bb[i,] + omega/sig2*eps
    vbb <- sb[i]^2 - (omega^2)/sig2
	vbb = ifelse(vbb<1*10^-32, .0001, vbb)
    s <- ifelse(vbb>=0 & vbb!=Inf & !is.na(vbb),sqrt(vbb),NaN)
    bounds <- bounds1(x,t,n)
    out<- NULL
    for(j in 1:length(x[ok])){
      out[ok][j] <- rtnorm(1, mean=mbb[ok][j], sd=s[ok][j],
                           lower=bounds[ok,][j,1],
                           upper=bounds[ok,][j,2])
    }
    out[wh] <- NA
    out[bl] <- t[bl]
    out[cT1] <- bounds[cT1,1]
    out[cT0] <- bounds[cT0,1]
    betab[,i] = out
  }
  omx <- 1 - x
  for (j in 1:length(x[ok])){
    betabs <- betab[ok,][j,]
    betaw[ok,][j,] <- t[ok][j]/omx[ok][j]-betabs*x[ok][j]/omx[ok][j]
  }

  if(sum(wh)>0){
    betaw[wh,] <- as.matrix(rep(1,dim(keep)[1]))%*%t(as.matrix(t[wh]))
  }

  if(sum(bl)>0){
    betaw[bl,] <- NA
    #betaw[bl,] <- as.matrix(rep(1,dim(keep)[1]))%*%t(as.matrix(t[bl]))
  }
  if(sum(cT1)>0){
    betaw[cT1,] <-
      as.matrix(rep(1,dim(keep)[1]))%*%t(as.matrix(bounds[cT1,3]))
  }
  if(sum(cT0)>0){
    betaw[cT0,]<-
      as.matrix(rep(1,dim(keep)[1]))%*%t(as.matrix(bounds[cT0,3]))
  }

  mbetab <- apply(betab,1,mean)
  mbetaw <- apply(betaw,1,mean)
  sdbetab <- apply(betab,1,sd)
  sdbetaw <- apply(betaw,1,sd)
  output <- list(ei.object$phi, ei.object$hessian, hessian, psi, mbetab, mbetaw,
                 sdbetab, sdbetaw, betab, betaw,resamp, erho, esigma,
                 ebeta, ealphab, ealphaw, numb, x, t, n, Zb, Zw,
                 truth,
                 precision, id)

  names(output) <- c("phi", "hessian", "hessianC","psi", "betab", "betaw",
                     "sbetab",
                     "sbetaw", "betabs", "betaws", "resamp", "erho",
                     "esigma", "ebeta", "ealphab", "ealphaw", "numb",
                     "x", "t", "n", "Zb", "Zw",
                     "truth", "precision", "id")
  class(output) <- "ei"
  return(output)
}




