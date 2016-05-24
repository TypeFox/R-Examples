#' Simulate closed population capture-mark-recapture data arising from multiple non-invasive marks
#'
#' This function generates encounter histories from simulated closed population capture-mark-recapture data consisting of multiple non-invasive marks. 
#'
#'
#' @param N True population size or abundance.
#' @param noccas The number of sampling occasions.
#' @param pbeta Logit- or probit-scale intercept term(s) for capture probability (p). Must be a scaler or vector of length \code{noccas}.
#' @param tau Additive logit- or probit-scale behavioral effect term for recapture probability (c).
#' @param sigma2_zp Logit- or probit-scale individual heterogeneity variance term.
#' @param delta_1 Conditional probability of type 1 encounter, given detection.
#' @param delta_2 Conditional probability of type 2 encounter, given detection.
#' @param alpha Conditional probability of simultaneous type 1 and type 2 detection, given both types encountered. Only applies when \code{data.type="sometimes"}.
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param link Link function for detection probability. Must be "\code{logit}" or "\code{probit}". Note that \code{\link{multimarkClosed}} is currently implemented for the logit link only.
#'
#' @return A list containing the following:
#' \item{Enc.Mat}{A matrix containing the observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions.}
#' \item{trueEnc.Mat}{A matrix containing the true (latent) encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions.}
#' @author Brett T. McClintock 
#' @seealso \code{\link{processdata}}, \code{\link{multimarkClosed}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' #simulate data for data.type="sometimes" using defaults
#' data<-simdataClosed(data.type="sometimes")
simdataClosed <- function(N=100,noccas=5,pbeta=-0.4,tau=0,sigma2_zp=0,delta_1=0.4,delta_2=0.4,alpha=0.5,data.type="never",link="logit"){
  
  if(length(pbeta)==1){
    pbeta=rep(pbeta,noccas)
  } else if(length(pbeta)!=noccas){
    stop(paste0("'pbeta' must be a scaler or vector of length ",noccas))
  }
  delta_B<-1-(delta_1+delta_2)
  if(delta_B<0) stop ("delta_1 and delta_2 must have sum less than 1")
  
  if(data.type=="never"){
    alpha<-0
  } else if(data.type=="always"){
    alpha<-1
  } else if(data.type!="sometimes"){
    stop("'data.type' must be 'never', 'sometimes', or 'always'")
  }
  
  tEnc.Mat<-matrix(0,nrow=N,ncol=noccas)        #"true" latent histories
  zp<-rnorm(N,0,sqrt(sigma2_zp))
  for(i in 1:N){
    ind<-0
    for(j in 1:noccas){
      if(link=="probit"){
        p<-pnorm(pbeta[j]+zp[i])
        c<-pnorm(pbeta[j]+tau+zp[i])
      } else if(link=="logit"){
        p<-expit(pbeta[j]+zp[i])
        c<-expit(pbeta[j]+tau+zp[i])        
      } else {stop("link function must be 'probit' or 'logit'")}
      tEnc.Mat[i,j] <- rbinom(1,1,((1-ind)*p+ind*c) )       #"true" latent histories
      if(tEnc.Mat[i,j]==1){
        ind<-1
      }
    }
  }
  Rand.Mat<-matrix(runif(N*noccas,0,1),N,noccas)
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat<delta_2)] <- 2      # type 2 encounters
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat>(1-delta_B))] <- 4  # type 1 and type 2 encounters
  tEnc.Mat[which(tEnc.Mat==4)] <- tEnc.Mat[which(tEnc.Mat==4)]-(runif(base::sum(tEnc.Mat==4))<(1-alpha))   # unobserved type 1 and type 2 encounters
  
  Enc.Mat <- get_Enc(tEnc.Mat,data.type)
  return(list(Enc.Mat=Enc.Mat,trueEnc.Mat=tEnc.Mat))
}

get_DMClosed<-function(mod.p,mod.delta,Enc.Mat,covs,type="Closed",...){
  Enc.Mat[which(Enc.Mat>0)] <- 1
  ch<-as.character(as.matrix( apply(Enc.Mat, 1, paste, collapse = ""), ncol=1 ))
  if(length(which(ch==paste(rep(0,ncol(Enc.Mat)),collapse="")))){
    ch<-ch[-which(ch==paste(rep(0,ncol(Enc.Mat)),collapse=""))] 
  }
  CH<-process.data(data.frame(ch,options(stringsAsFactors = FALSE)),model=type)
  temp<-make.design.data(CH,...)
  if(CH$nocc %% 2 == 0){
    temp$p$Time <- temp$p$Time+1-(CH$nocc/2)-.5
  } else {
    temp$p$Time <- temp$p$Time+1-ceiling(CH$nocc/2)
  }
  temp$c<-temp$p
  temp$c$model.index<-temp$c$model.index+CH$nocc
  temp$c$c<-1
  modterms<-attributes(terms(mod.p))$term.labels
  if(any(modterms=="time")){
    mod.p<-formula(paste("~-1",paste(modterms,collapse="+"),sep="+"))
  }
  if(any(modterms=="h")){
    if(length(modterms)>1){
      mod.p<-formula(paste("~-1",paste(modterms[-which(modterms=="h")],collapse="+"),sep="+"))
    } else {
      mod.p<-formula(~1)  
    }
    mod.p.h<-TRUE 
  } else {
    mod.p.h<-FALSE
  }
  if(length(covs)){
    temp$p<-cbind(temp$p,covs)
    temp$c<-cbind(temp$c,covs)
  }
  DMp<-model.matrix(mod.p,temp$p)
  DMc<-model.matrix(mod.p,temp$c)
  
  if(nrow(DMp)!=CH$nocc) stop(paste("model design matrix must have",CH$nocc,"rows"))
  
  if(!any(colnames(DMp)=="(Intercept)")){
    DMp<-cbind(rep(1,CH$nocc),DMp)
    colnames(DMp)[1]<-"(Intercept)"
    DMc<-cbind(rep(1,CH$nocc),DMc)
    colnames(DMc)[1]<-"(Intercept)"
  }
  rownames(DMp) <- paste0("p[",1:CH$nocc,"]")
  rownames(DMc) <- paste0("c[",1:CH$nocc,"]")
  
  deltattr <- attributes(terms(mod.delta))$term.labels
  if(length(deltattr)){
    if(deltattr!="type") stop("'mod.delta' must be '~1' or '~type'")
  }
  return(list(p=DMp,c=DMc,mod.p.h=mod.p.h,mod.delta=mod.delta))
}

pstarintegrand<-function(beta,sigma,DM,gq){
  
  nodes <- gq$nodes
  weights <- gq$weights
  npts <- length(nodes)
  
  XB <- DM %*% beta
  temp = apply(matrix(1. - expit(sqrt(2.0)*sigma*rep(nodes,nrow(DM))+rep(XB,each=npts)),nrow=npts),1,prod)
  oneminuspstar <- sum(1./sqrt(pi)*weights*temp)
  1.-min(1.-tol,max(tol,oneminuspstar))
}

mcmcClosed<-function(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,gq,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog){
  
  weights<-gq$weights
  nodes<-gq$nodes
  npoints<-length(weights)
  
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  DMp<-DM$p
  DMc<-DM$c
  mod.p.h<-DM$mod.p.h
  pdim<-ncol(DMp)
  
  #declare and initialize parameters
  pbeta<-rep(NA,max(1,floor(iter/thin))*(pdim))
  zp<-rep(NA,ifelse(any(params=="zp"),max(1,floor(iter/thin))*M,M))
  H<-rep(NA,ifelse(any(params=="H"),max(1,floor(iter/thin))*M,M))
  sigma2_zp<-rep(NA,max(1,floor(iter/thin)))
  alpha<-rep(NA,max(1,floor(iter/thin)))
  delta_1<-rep(NA,max(1,floor(iter/thin)))
  delta_2<-rep(NA,max(1,floor(iter/thin)))
  N<-rep(NA,max(1,floor(iter/thin)))
  psi<-rep(NA,max(1,floor(iter/thin)))
  logPosterior<-rep(NA,max(1,floor(iter/thin)))
  
  pbeta[1:pdim] <- inits[[ichain]]$pbeta
  zp[1:M] <- inits[[ichain]]$zp
  H[1:M] <- inits[[ichain]]$H-1
  sigma2_zp[1] <- inits[[ichain]]$sigma2_zp
  alpha[1] <- inits[[ichain]]$alpha
  delta_1[1] <- inits[[ichain]]$delta_1
  delta_2[1] <- inits[[ichain]]$delta_2
  N[1] <- inits[[ichain]]$N
  psi[1] <- inits[[ichain]]$psi
  
  arate<-numeric(M+pdim+1)
  
  posterior <- .C(ClosedC,as.integer(ichain),as.numeric(mu0), as.numeric(sigma2_mu0), as.numeric(pbeta), as.numeric(zp), as.numeric(sigma2_zp), as.numeric(delta_1),as.numeric(delta_2),as.numeric(alpha), as.integer(inits[[ichain]]$x), as.numeric(N), as.numeric(psi), as.integer(H),
                  as.integer(noccas), as.integer(M), as.numeric(a0delta), as.numeric(a0alpha), as.numeric(b0alpha), as.numeric(a), as.numeric(a0psi), as.numeric(b0psi),
                  as.numeric(Prop.sd),as.numeric(arate),as.numeric(logPosterior),
                  as.integer(length(mms@vAll.hists)/noccas),as.integer(mms@vAll.hists), as.integer(mms@C), as.integer(mms@indBasis-1), as.integer(mms@ncolbasis), as.integer(mms@knownx), as.numeric(as.vector(t(DMp))), as.numeric(as.vector(t(DMc))),as.integer(pdim),
                  as.integer(iter), as.integer(thin), as.integer(adapt), as.integer(bin), as.numeric(taccept),as.numeric(tuneadjust),as.integer(maxnumbasis),
                  as.integer(npoints),as.numeric(weights),as.numeric(nodes),as.integer(mod.p.h),as.integer(mms@data.type=="sometimes"),as.integer(any(params=="zp")),as.integer(any(params=="H")),as.integer(DM$mod.delta != ~NULL),as.integer(DM$mod.delta==formula(~type)),as.integer(printlog),NAOK = TRUE) 
  
  names(posterior) <- c("ichain","mu_0","sigma2_mu","pbeta", "zp", "sigma2_zp", "delta_1","delta_2","alpha", "x", "N", "psi","H", "noccas", "M","a0delta", "a0alpha", "b0alpha","a","a0psi","b0psi","Prop.sd","arate","logPosterior","nHists","vAll.hists","C", "indBasis", "ncolBasis","knownx","DMp","DMc","pdim","iter", "thin", "adapt", "bin", "taccept","tuneadjust","maxnumbasis","npoints","weights","nodes","mod.p.h","sometimes?","zp?","H?","updatedelta?","type?","printlog?")
  
  g <- posterior$iter
  x <- posterior$x
  if(any(params=="zp")){
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),matrix(posterior$zp[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)],ncol=M,byrow=T),posterior$sigma2_zp[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    zp <- NULL
  } else {
    zp <- posterior$zp
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),posterior$sigma2_zp[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  if(any(params=="H")){
    posterior<-cbind(temp,matrix(posterior$H[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)]+1,ncol=M,byrow=T),posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    H <- NULL
  } else {
    H <- posterior$H+1
    posterior<-cbind(temp,posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  return(list(posterior=posterior,x=x,H=H,zp=zp,g=g))
}

loglikeClosed<-function(parms,DM,noccas,C,All.hists,gq){
  
  H <- parms$H
  pbeta <- parms$pbeta
  if(DM$mod.p.h){
    zp <- parms$zp
  } else {
    zp <- rep(0,length(H))
  }
  if(DM$mod.delta != ~NULL){
    if(DM$mod.delta==formula(~type)){
      delta_1 <- parms$delta_1
      delta_2 <- parms$delta_2
    } else {
      delta_1 <- delta_2 <- parms$delta
    }
    alpha <- parms$alpha
  } else {
    delta_1 <- 1.0
    delta_2 <- 0.0
    alpha <- 0.0
  }
  
  Hind <- H[which(H>1)]
  indhist <- All.hists[Hind,]
  n<-length(Hind)
  firstcap<- (C[Hind]>=matrix(rep(1:noccas,each=n),nrow=n,ncol=noccas))
  
  p <- expittol(matrix(rep(DM$p%*%pbeta,each=n)*firstcap+rep(DM$c%*%pbeta,each=n)*(1-firstcap)+zp[which(H>1)],nrow=n,ncol=noccas))
  
  loglike <- base::sum( log( (indhist==0) * (1. - p)
                             + (indhist==1) * p * delta_1  
                             + (indhist==2) * p * delta_2
                             + (indhist==3) * p * (1. - delta_1 - delta_2) * (1. - alpha)
                             + (indhist==4) * p * (1. - delta_1 - delta_2) * alpha ))
  if(DM$mod.p.h){
    pstar <- pstarintegrand(pbeta,sqrt(parms$sigma2_zp),DM$p,gq)
  } else {
    pstar <- 1-min(1.-tol,max(tol,prod(1-expit(DM$p%*%pbeta))))
  }    
  loglike <- loglike + dbinom(n,parms$N,pstar,1) - n * log(pstar)  
  loglike
}

priorsClosed<-function(parms,DM,priorparms,data_type){
  
  priors <- (base::sum(dnorm(parms$pbeta,priorparms$mu0,sqrt(priorparms$sigma2_mu0),log=TRUE))
             + -log(parms$N))
  
  if(DM$mod.delta != ~NULL){
    if(DM$mod.delta==formula(~type)){
      priors <- priors + ddirichlet(c(parms$delta_1,parms$delta_2,1.-parms$delta_1-parms$delta_2),priorparms$a0delta)
    } else {
      priors <- priors + dbeta(2*parms$delta,priorparms$a0delta[1],priorparms$a0delta[2],log=TRUE)
    }
    if(data_type=="sometimes"){
      priors <- priors + dbeta(parms$alpha,priorparms$a0alpha,priorparms$b0alpha,log=TRUE)
    }
    priors <- priors + (base::sum(dbinom((parms$H>1),1,parms$psi,log=TRUE))
                         + dbeta(parms$psi,priorparms$a0psi,priorparms$b0psi,log=TRUE))
  }
  
  if(DM$mod.p.h){
    priors <- priors + (base::sum(dnorm(parms$zp,0.0,sqrt(parms$sigma2_zp),log=TRUE))
                        + log(2.0*dcauchy(sqrt(parms$sigma2_zp),0.0,priorparms$a,log=FALSE)))
  }        
  priors
}

posteriorClosed<-function(parms,DM,mms,priorparms,gq){
  nchains<-length(parms)
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas)
  for(ichain in 1:nchains){
    temp<-parms[[ichain]]
    
    loglike <- loglikeClosed(temp,DM,noccas,mms@C,All.hists,gq)
    
    if(!is.finite(loglike)) {
      stop(paste0("initial model likelihood is ",loglike," for chain ",ichain,". Try different initial values."))
    }
    
    posterior <- loglike + priorsClosed(temp,DM,priorparms,mms@data.type)
    
    if(!is.finite(posterior)) {
      stop(paste("initial model posterior is",posterior,"for chain",ichain,". Try different initial values or prior parameters"))
    }
  }
}

checkClosed<-function(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,npoints,maxnumbasis,a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi){
  
  if(mms@data.type!="sometimes" & any(parms=="alpha")) stop("Parameter 'alpha' only applies to models for the 'sometimes' data type")
  
  params<-parms
  if(any(parms=="all")){
    if(mms@data.type=="sometimes"){
      params<-parmlist
    } else {
      params<-parmlist[which(parmlist!="alpha")]
    }
  } else {
    if(!all(match(params,parmlist,nomatch=0))) stop(paste0("'",params[match(params,parmlist,nomatch=0)==0],"' is not a valid parameter\n  "))
  }  
  
  if(adapt<0) stop("'adapt' must be >=0")
  if((bin<1 | bin>iter) & iter>0) stop("'bin' must be >0 and <",iter)
  if(thin>max(1,floor((iter-burnin+1)/2)) | thin<1) stop("'thin' must be >0 and <=",max(1,floor((iter-burnin+1)/2)))
  if(taccept<=0 | taccept>1) stop ("'taccept' must be >0 and <=1")
  if(tuneadjust<=0 | tuneadjust>1) stop ("'tuneadjust' must be >0 and <=1")
  if(npoints<1) stop("'npoints' must be greater than 0")
  if(mms@ncolbasis & (maxnumbasis<1 | maxnumbasis>mms@ncolbasis)) stop("'maxnumbasis' must be between 1 and ",mms@ncolbasis)
  if(!all(c(a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi)>0)) stop("'a0delta', 'a0alpha', 'b0alpha', 'a', 'sigma2_mu0', 'a0psi', and 'b0psi' must be >0")
  
  mod.p.h<-DM$mod.p.h
  if(any(parms=="all")){
    if(!mod.p.h){
      params<-params[which(params!="zp" & params!="sigma2_zp")]
    }
  } else {
    if(!mod.p.h & (any(params=="zp") | any(params=="sigma2_zp"))) stop("Parameters 'sigma2_zp' and 'zp' only apply to individual heterogeneity models")
  }
  pdim<-ncol(DM$p)
  if(!pdim) stop("'mod.p' must include at least 1 parameter")
  
  params
}

processClosedchains<-function(chains,params,DM,M,noccas,nchains,iter,burnin,thin){
  
  parms<-params
  if(any(parms=="pbeta")){
    parms<-c(paste0("pbeta[",colnames(DM$p),"]"),params[which(params!="pbeta")])
  }
  if(any(parms=="delta")){
    if(DM$mod.delta==formula(~type)){
      deltaname<-c("delta_1","delta_2")   
      parms<-c(parms[which(parms!="delta")],"delta_1","delta_2") 
    } else {
      deltaname<-c("delta")   
      parms<-c(parms[which(parms!="delta")],"delta_1") 
    }
  } else {
    deltaname<-NULL
  }
  if(any(parms=="zp")){
    zpname<-paste0("zp[",1:M,"]")
    parms<-c(zpname,parms[which(parms!="zp")])
  } else {
    zpname<-NULL
  }
  if(any(parms=="H")){
    Hname<-paste0("H[",1:M,"]")
    parms<-c(Hname,parms[which(parms!="H")])
  } else {
    Hname<-NULL
  }
  
  initial.values <- list()
  
  for(ichain in 1:nchains){
    checkend <- chains[[ichain]]$g
    if(checkend<iter | !is.finite(chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),ncol(chains[[ichain]]$posterior)])) {
      warning(paste0("chain ",ichain," terminated at iteration ",checkend,"; check log for more information"))
      if(!checkend & burnin<1){
        initstemp <- chains[[ichain]]$posterior[1,]
      } else if(floor(checkend/thin)>floor(burnin/thin)){
        initstemp <- chains[[ichain]]$posterior[floor(checkend/thin)-floor(burnin/thin),]  
      } else {
        initstemp <- chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),]
      }
    } else {
      initstemp <- chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),]
    }
    names(initstemp) <- c(paste0("pbeta[",colnames(DM$p),"]"),zpname,"sigma2_zp","delta_1","delta_2","alpha","N","psi",Hname,"logPosterior")
    if(any(params=="zp")){
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],zp=initstemp[zpname],sigma2_zp=initstemp["sigma2_zp"],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],N=initstemp["N"],psi=initstemp["psi"],x=chains[[ichain]]$x,H=chains[[ichain]]$H)
    } else {
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],zp=chains[[ichain]]$zp,sigma2_zp=initstemp["sigma2_zp"],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],N=initstemp["N"],psi=initstemp["psi"],x=chains[[ichain]]$x,H=chains[[ichain]]$H)
      names(initial.values[[ichain]]$zp) <- paste0("zp[",1:M,"]")
    }
    if(any(params=="H")){
      initial.values[[ichain]]$H <- initstemp[Hname]
    } else {
      initial.values[[ichain]]$H <- chains[[ichain]]$H
      names(initial.values[[ichain]]$H) <- paste0("H[",1:M,"]")
    }
    names(initial.values[[ichain]]$x) <- paste0("x[",1:length(initial.values[[ichain]]$x),"]")
    chains[[ichain]] <- chains[[ichain]]$posterior
    colnames(chains[[ichain]]) <- names(initstemp)   
    chains[[ichain]] <- chains[[ichain]][,parms]
    if(!is.null(deltaname)){
      if(!is.null(nrow(chains[[ichain]]))) {
        colnames(chains[[ichain]])[which(substr(colnames(chains[[ichain]]),1,nchar("delta"))=="delta")] <- deltaname
      } else {
        names(chains[[ichain]])[which(substr(names(chains[[ichain]]),1,nchar("delta"))=="delta")] <- deltaname     
      }
    }
    chains[[ichain]] <- mcmc(chains[[ichain]],start=1,thin=1)
    if(burnin<thin){
      temp=seq(thin,max(1,iter),thin)
    } else {
      temp=seq(thin*(floor(burnin/thin)+1),iter,thin)
    }
    attributes(chains[[ichain]])$mcpar <- c(head(temp,n=1),tail(temp,n=1),thin)  
  }
  chains<-as.mcmc.list(chains)
  return(list(chains=chains,initial.values=initial.values))  
}

#' Fit closed population abundance models for ``traditional'' capture-mark-recapture data consisting of a single mark type
#'
#' This function fits closed population abundance models for ``traditional'' capture-mark-recapture data consisting of a single mark type using Bayesian analysis methods. Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix of observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions. With a single mark type, encounter histories consist of only non-detections (0) and type 1 encounters (1).
#' @param covs A data frame of temporal covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of sampling occasions. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mod.p Model formula for detection probability. For example, \code{mod.p=~1} specifies no effects (i.e., intercept only), \code{mod.p~time} specifies temporal effects, \code{mod.p~c} specifies behavioral reponse (i.e., trap "happy" or "shy"), \code{mod.p~h} specifies individual heterogeneity, and \code{mod.p~time+c} specifies additive temporal and behavioral effects.
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are logit-scale detection probability parameters ("\code{pbeta}"), population abundance ("\code{N}"), logit-scale individual heterogeneity variance term ("\code{sigma2_zp}"), and logit-scale individual effects ("\code{zp}"). The log posterior density ("\code{logPosterior}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt The number of iterations for proposal distribution adaptation. If \code{adapt = 0} then no adaptation occurs.
#' @param bin Bin length for calculating acceptance rates during adaptive phase (\code{0 < bin <= iter}).
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Target acceptance rate during adaptive phase (\code{0 < taccept <= 1}). Acceptance rate is monitored every \code{bin} iterations. Default is \code{taccept = 0.44}.
#' @param tuneadjust Adjustment term during adaptive phase (\code{0 < tuneadjust <= 1}). If acceptance rate is less than \code{taccept}, then proposal term (\code{proppbeta}, \code{propzp}, or \code{propsigmap}) is multiplied by \code{tuneadjust}. If acceptance rate is greater than or equal to \code{taccept}, then proposal term is divided by \code{tuneadjust}. Default is \code{tuneadjust = 0.95}.
#' @param proppbeta Scaler or vector (of length k) specifying the initial standard deviation of the Normal(pbeta[j], proppbeta[j]) proposal distribution. If \code{proppbeta} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{proppbeta = 0.1}.
#' @param propzp Scaler or vector (of length M) specifying the initial standard deviation of the Normal(zp[i], propzp[i]) proposal distribution. If \code{propzp} is a scaler, then this value is used for all i = 1, ..., M individuals. Default is \code{propzp = 1}.
#' @param propsigmap Scaler specifying the initial Gamma(shape = 1/\code{propsigmap}, scale = sigma_zp * \code{propsigmap}) proposal distribution for sigma_zp = sqrt(sigma2_zp). Default is \code{propsigmap=1}.
#' @param npoints Number of Gauss-Hermite quadrature points to use for numerical integration. Accuracy increases with number of points, but so does computation time.
#' @param a Scale parameter for [sigma_z] ~ half-Cauchy(a) prior for the individual hetegeneity term sigma_zp = sqrt(sigma2_zp). Default is ``uninformative'' \code{a = 25}.
#' @param mu0 Scaler or vector (of length k) specifying mean of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{mu0 = 0}.
#' @param sigma2_mu0 Scaler or vector (of length k) specifying variance of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{sigma2_mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{sigma2_mu0 = 1.75}.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for "\code{pbeta}", "\code{zp}", "\code{sigma2_zp}", and "\code{N}". Default is \code{initial.values = NULL}, which causes initial values to be generated automatically.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p}. See \code{\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{markClosed} (or \code{\link{markCJS}}) is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{markClosed}. Note that setting \code{parms="all"} is required for any \code{markClosed} model output to be used in \code{\link{multimodelClosed}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.delta}{Formula always \code{NULL}; only for internal use in \code{\link{multimodelClosed}}.}
#' \item{DM}{A list of design matrices for detection probability generated for model \code{mod.p}, where DM$p is the design matrix for initial capture probability (p) and DM$c is the design matrix for recapture probability (c).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{zp}", "\code{sigma2_zp}", and "\code{N}".}
#' \item{mms}{An object of class \code{multimarksetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimodelClosed}}
#' @examples
#' \dontshow{
#' test<-markClosed(Enc.Mat=simdataClosed(delta_1=1,delta_2=0)$Enc.Mat,iter=10,burnin=0,bin=5)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run single chain using the default model for simulated ``traditional'' data
#' data<-simdataClosed(delta_1=1,delta_2=0)$Enc.Mat
#' sim.dot<-markClosed(data)
#' 
#' #Posterior summary for monitored parameters
#' summary(sim.dot$mcmc)
#' plot(sim.dot$mcmc)}
markClosed<-function(Enc.Mat,covs=data.frame(),mod.p=~1,parms=c("pbeta","N"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propzp=1,propsigmap=1,npoints=500,a=25,mu0=0,sigma2_mu0=1.75,initial.values=NULL,printlog=FALSE,...){
  if(any(Enc.Mat>1 | Enc.Mat<0)) stop("With a single mark type, encounter histories can only contain 0's (non-detections) and 1's (detections)")
  mms <- processdata(Enc.Mat,covs=covs,known=rep(1,nrow(Enc.Mat)))
  out <- multimarkClosed(mms=mms,mod.p=mod.p,mod.delta=~NULL,parms=parms,nchains=nchains,iter=iter,adapt=adapt,bin=bin,thin=thin,burnin=burnin,taccept=taccept,tuneadjust=tuneadjust,proppbeta=proppbeta,propzp=propzp,propsigmap=propsigmap,npoints=npoints,a=a,mu0=mu0,sigma2_mu0=sigma2_mu0,initial.values=initial.values,printlog=printlog,...)
  out$initial.values <- lapply(out$initial.values,function(x) list(pbeta=x$pbeta,zp=x$zp,sigma2_zp=x$sigma2_zp,N=x$N))
  return(out)
}

#' Fit closed population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks
#'
#' This function fits closed population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks using Bayesian analysis methods. Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix of observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions (ignored unless \code{mms=NULL}).
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param covs A data frame of temporal covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of sampling occasions. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mms An optional object of class \code{multimarksetup-class}; if \code{NULL} it is created. See \code{\link{processdata}}.
#' @param mod.p Model formula for detection probability. For example, \code{mod.p=~1} specifies no effects (i.e., intercept only), \code{mod.p~time} specifies temporal effects, \code{mod.p~c} specifies behavioral reponse (i.e., trap "happy" or "shy"), \code{mod.p~h} specifies individual heterogeneity, and \code{mod.p~time+c} specifies additive temporal and behavioral effects.
#' @param mod.delta Model formula for conditional probabilities of type 1 (delta_1) and type 2 (delta_2) encounters, given detection. Currently only \code{mod.delta=~1} (i.e., \eqn{\delta_1 = \delta_2}) and \code{mod.delta=~type} (i.e., \eqn{\delta_1 \ne \delta_2}) are implemented.
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are logit-scale detection probability parameters ("\code{pbeta}"), population abundance ("\code{N}"), conditional probability of type 1 or type 2 encounter, given detection ("\code{delta})", probability of simultaneous type 1 and type 2 detection, given both types encountered ("\code{alpha}"), logit-scale individual heterogeneity variance term ("\code{sigma2_zp}"), logit-scale individual effects ("\code{zp}"), and the probability that a randomly selected individual from the \code{M = nrow(Enc.Mat)} observed individuals belongs to the \eqn{n} unique individuals encountered at least once ("\code{psi}"). Individual encounter history indices ("\code{H}") and the log posterior density ("\code{logPosterior}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt The number of iterations for proposal distribution adaptation. If \code{adapt = 0} then no adaptation occurs.
#' @param bin Bin length for calculating acceptance rates during adaptive phase (\code{0 < bin <= iter}).
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Target acceptance rate during adaptive phase (\code{0 < taccept <= 1}). Acceptance rate is monitored every \code{bin} iterations. Default is \code{taccept = 0.44}.
#' @param tuneadjust Adjustment term during adaptive phase (\code{0 < tuneadjust <= 1}). If acceptance rate is less than \code{taccept}, then proposal term (\code{proppbeta}, \code{propzp}, or \code{propsigmap}) is multiplied by \code{tuneadjust}. If acceptance rate is greater than or equal to \code{taccept}, then proposal term is divided by \code{tuneadjust}. Default is \code{tuneadjust = 0.95}.
#' @param proppbeta Scaler or vector (of length k) specifying the initial standard deviation of the Normal(pbeta[j], proppbeta[j]) proposal distribution. If \code{proppbeta} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{proppbeta = 0.1}.
#' @param propzp Scaler or vector (of length M) specifying the initial standard deviation of the Normal(zp[i], propzp[i]) proposal distribution. If \code{propzp} is a scaler, then this value is used for all i = 1, ..., M individuals. Default is \code{propzp = 1}.
#' @param propsigmap Scaler specifying the initial Gamma(shape = 1/\code{propsigmap}, scale = sigma_zp * \code{propsigmap}) proposal distribution for sigma_zp = sqrt(sigma2_zp). Default is \code{propsigmap=1}.
#' @param npoints Number of Gauss-Hermite quadrature points to use for numerical integration. Accuracy increases with number of points, but so does computation time.
#' @param maxnumbasis Maximum number of basis vectors to use when proposing latent history frequency updates. Default is \code{maxnumbasis = 1}, but higher values can potentially improve mixing.
#' @param a0delta Scaler or vector (of length d) specifying the prior for the conditional (on detection) probability of type 1 (delta_1), type 2 (delta_2), and both type 1 and type 2 encounters (1-delta_1-delta_2). If \code{a0delta} is a scaler, then this value is used for all a0delta[j] for j = 1, ..., d. For \code{mod.delta=~type}, d=3 with [delta_1, delta_2, 1-delta_1-delta_2] ~ Dirichlet(a0delta) prior. For \code{mod.delta=~1}, d=2 with [tau] ~ Beta(a0delta[1],a0delta[2]) prior, where (delta_1,delta_2,1-delta_1-delta_2) = (tau/2,tau/2,1-tau). See McClintock et al. (2013) for more details.
#' @param a0alpha Specifies "shape1" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{a0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param b0alpha Specifies "shape2" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{b0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param a Scale parameter for [sigma_z] ~ half-Cauchy(a) prior for the individual hetegeneity term sigma_zp = sqrt(sigma2_zp). Default is ``uninformative'' \code{a = 25}.
#' @param mu0 Scaler or vector (of length k) specifying mean of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{mu0 = 0}.
#' @param sigma2_mu0 Scaler or vector (of length k) specifying variance of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{sigma2_mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{sigma2_mu0 = 1.75}.
#' @param a0psi Specifies "shape1" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{a0psi = 1}.
#' @param b0psi Specifies "shape2" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{b0psi = 1}.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for parameters and latent variables. Default is \code{initial.values = NULL}, which causes initial values to be generated automatically. In addition to the parameters ("\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_zp}", "\code{zp}", and "\code{psi}"), initial values can be specified for the initial latent history frequencies ("\code{x}") and initial individual encounter history indices ("\code{H}").
#' @param known Optional integer vector indicating whether the encounter history of an individual is known with certainty (i.e., the observed encounter history is the true encounter history). Encounter histories with at least one type 4 encounter are automatically assumed to be known, and \code{known} does not need to be specified unless there exist encounter histories that do not contain a type 4 encounter that happen to be known with certainty (e.g., from independent telemetry studies). If specified, \code{known = c(v_1,v_2,...,v_M)} must be a vector of length \code{M = nrow(Enc.Mat)} where \code{v_i = 1} if the encounter history for individual \code{i} is known (\code{v_i = 0} otherwise). Note that known all-zero encounter histories (e.g., `000') are ignored.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p}. See \code{\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{multimarkClosed} (or \code{\link{multimarkCJS}}) is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{multimarkClosed}. Note that setting \code{parms="all"} is required for any \code{multimarkClosed} model output to be used in \code{\link{multimodelClosed}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.delta}{Model formula for conditional probability of type 1 or type 2 encounter, given detection (as specified by \code{mod.delta} above).}
#' \item{DM}{A list of design matrices for detection probability generated for model \code{mod.p}, where DM$p is the design matrix for initial capture probability (p) and DM$c is the design matrix for recapture probability (c).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_zp}", "\code{zp}", "\code{psi}", "\code{x}", and "\code{H}".}
#' \item{mms}{An object of class \code{multimarksetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{bobcat}}, \code{\link{processdata}}, \code{\link{multimodelClosed}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' 
#' McClintock, B. T., Bailey, L. L., Dreher, B. P., and Link, W. A. 2014. Probit models for capture-recapture data subject to imperfect detection, individual heterogeneity and misidentification. \emph{The Annals of Applied Statistics} 8: 2461-2484.
#' @examples
#' \dontshow{
#' test<-multimarkClosed(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0,bin=5)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run single chain using the default model for bobcat data
#' bobcat.dot<-multimarkClosed(bobcat)
#' 
#' #Posterior summary for monitored parameters
#' summary(bobcat.dot$mcmc)
#' plot(bobcat.dot$mcmc)}
multimarkClosed<-function(Enc.Mat,data.type="never",covs=data.frame(),mms=NULL,mod.p=~1,mod.delta=~type,parms=c("pbeta","delta","N"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propzp=1,propsigmap=1,npoints=500,maxnumbasis=1,a0delta=1,a0alpha=1,b0alpha=1,a=25,mu0=0,sigma2_mu0=1.75,a0psi=1,b0psi=1,initial.values=NULL,known=integer(),printlog=FALSE,...){
  
  if(is.null(mms)) mms <- processdata(Enc.Mat,data.type,covs,known)
  if(class(mms)!="multimarksetup") stop("'mms' must be an object of class 'multimarksetup'")
  validObject(mms)
  
  if(class(mod.p)!="formula") stop("'mod.p' must be an object of class 'formula'")
  if(class(mod.delta)!="formula") stop("'mod.delta' must be an object of class 'formula'")
  DM<-get_DMClosed(mod.p,mod.delta,mms@Enc.Mat,covs=mms@covs,...)
  
  if(iter>0){
    if(iter<=burnin) stop(paste("'burnin' must be less than ",iter))
  } else {
    burnin<-0
  }
  
  if(mod.delta != ~NULL) {
    parmlist<-c("pbeta","delta","N","sigma2_zp","zp","alpha","psi","H","logPosterior")
  } else {
    parmlist<-c("pbeta","N","sigma2_zp","zp","logPosterior")    
  }
  params <- checkClosed(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,npoints,maxnumbasis,a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi)
  
  data.type<-mms@data.type
  Enc.Mat<-mms@Enc.Mat
  M<-nrow(Enc.Mat)
  noccas<-ncol(Enc.Mat)
  covs<-mms@covs
  pdim<-ncol(DM$p)
  
  mu0 <- checkvecs(mu0,pdim,"mu0")
  sigma2_mu0 <- checkvecs(sigma2_mu0,pdim,"sigma2_mu0")
  a0delta <- checkvecs(a0delta,ifelse(mod.delta==formula(~type),3,2),"a0delta")
  
  gq<-gauss.quad(npoints,kind="hermite")
  
  inits<-get_inits(mms,nchains,initial.values,M,data.type,a0alpha,b0alpha,a0delta,a0psi,b0psi,DM,gq)
  
  priorparms <-list(a0delta=a0delta,a0alpha=a0alpha,b0alpha=b0alpha,a=a,mu0=mu0,sigma2_mu0=sigma2_mu0,a0psi=a0psi,b0psi=b0psi,npoints=npoints)
  
  message("\nFitting closed capture model with logit link\n")
  if(mod.delta != ~NULL) message("data type = \"",data.type,"\"\n")
  message("p model = ",as.character(mod.p))
  if(mod.delta != ~NULL) message("delta model = ",as.character(mod.delta))
  message("\nInitializing model \n")
  posteriorClosed(inits,DM,mms,priorparms,gq)
  
  propzp <- checkvecs(propzp,M,"propzp")
  proppbeta <- checkvecs(proppbeta,pdim,"proppbeta")
  if(length(propsigmap)!=1) stop("'propsigmap' must be a scaler")
  
  Prop.sd <- c(propzp,proppbeta,propsigmap)
  
  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    modlog <- ifelse(mod.delta != ~NULL,"multimarkClosed","markClosed")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0(modlog,"_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("mcmcClosed"),envir=environment())  
    chains <- parLapply(cl,1:nchains, function(ichain) mcmcClosed(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,gq,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog))
    stopCluster(cl)
    gc()
  } else {
    chains <- vector('list',nchains)
    chains[[nchains]] <- mcmcClosed(nchains,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,gq,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog)
    gc()
  }
  
  chains <- processClosedchains(chains,params,DM,M,noccas,nchains,iter,burnin,thin)
  return(list(mcmc=chains$chains,mod.p=mod.p,mod.delta=mod.delta,DM=list(p=DM$p,c=DM$c),initial.values=chains$initial.values,priorparms=priorparms,mms=mms))
}

#' Calculate posterior capture and recapture probabilities
#'
#' This function calculates posterior capture (\eqn{p}) and recapture (\eqn{c}) probabilities for each sampling occasion from \code{\link{multimarkClosed}} output. 
#'
#'
#' @param out List of output returned by \code{\link{multimarkClosed}}.
#' @param link Link function for detection probability. Must be "\code{logit}" or "\code{probit}". Note that \code{\link{multimarkClosed}} is currently implemented for the logit link only.
#' @return An object of class \code{\link[coda]{mcmc.list}} containing the following:
#' \item{p}{Posterior samples for capture probability (\eqn{p}) for each sampling occasion.}
#' \item{c}{Posterior samples for recapture probability (\eqn{c}) for each sampling occasion.}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkClosed}}
#' @examples
#' \dontshow{
#' test<-getprobsClosed(multimarkClosed(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0,bin=5))}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run behavior model for bobcat data with constant detection probability (i.e., mod.p=~c)
#' bobcat.c <- multimarkClosed(bobcat,mod.p=~c)
#'   
#' #Calculate capture and recapture probabilities
#' pc <- getprobsClosed(bobcat.c)
#' summary(pc)}
getprobsClosed<-function(out,link="logit"){
  
  DMp<-out$DM$p
  DMc<-out$DM$c
  mod.p.h<-any("h"==attributes(terms(out$mod.p))$term.labels)
  
  noccas<-nrow(DMp)
  if(noccas<2) stop("must have >1 sampling occasion")
  
  pbetanames<-paste0("pbeta[",colnames(DMp),"]")
  nchains<-length(out$mcmc)
  
  pc<-vector("list",nchains)
  
  varind <- is.null(varnames(out$mcmc))
  if(!varind){
    vars <- varnames(out$mcmc)
  } else {
    vars <- names(out$mcmc[[1]])    
  }
  if(!any(match(pbetanames,vars,nomatch=0))) stop("'pbeta' parameters not found")
  
  for(ichain in 1:nchains){
    
    p <- inverseXB(ichain,out,pbetanames,mod.p.h,DMp,noccas,varind,vars,"p","sigma2_zp",link)
    rc <- inverseXB(ichain,out,pbetanames,mod.p.h,DMc,noccas,varind,vars,"c","sigma2_zp",link)
    
    if(dim(rc)[1]==1){
      rc <- matrix(rc[,-1],nrow=1)
    } else if(noccas<3){
      rc <- matrix(rc[,-1],ncol=1)      
    } else {
      rc <- rc[,-1]
    }
    colnames(p) <- paste0("p[",1:noccas,"]")
    colnames(rc) <- paste0("c[",2:noccas,"]")
    pc[[ichain]]<- mcmc(cbind(p,rc),start=start(out$mcmc),end=end(out$mcmc),thin=attributes(out$mcmc[[ichain]])$mcpar[3])
  }
  return(as.mcmc.list(pc))
}

checkparmsClosed <- function(mms,modlist,params,parmlist,M){    
  deltatypeind <- which(lapply(modlist,function(x) any("~type"==x$mod.delta))==TRUE)
  if(length(deltatypeind)){
    if(!all(lapply(params[deltatypeind],function(x) base::sum(match(x,c("delta_1","delta_2"),nomatch=0)))==base::sum(1:length(1:2)))) stop("required parameters not found for all models")
  }
  delta1ind <- which(lapply(modlist,function(x) any("~1"==x$mod.delta))==TRUE)
  if(length(delta1ind)){
    if(!all(lapply(params[delta1ind],function(x) base::sum(match(x,"delta",nomatch=0)))==1)) stop("required parameters not found for all models")
  }
  if(length(deltatypeind) | length(delta1ind)){
    parmlist<-c(parmlist,"psi",paste0("H[",1:M,"]"))
    if(mms@data.type=="sometimes"){
      parmlist<-c(parmlist,"alpha")
    }
    if((length(deltatypeind)+length(delta1ind))!=length(modlist)) stop("Cannot perform multimodel inference using both 'multimarkClosed()' and 'markClosed()' models")
  }
  hind <- which(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.p))$term.labels))==TRUE)  
  if(!length(hind)){
    if(!all(lapply(params,function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
  } else {
    if(!all(lapply(params[-hind],function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
    parmlist<-c(parmlist,"sigma2_zp",paste0("zp[",1:M,"]"))
    if(!all(lapply(params[hind],function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist))))  stop("required parameters not found for all models")
  }
}

monitorparmsClosed <- function(parms,parmlist,noccas){ 
  
  if(!all(match(parms,parmlist,nomatch=0))) stop(paste0("monitored parameters ('monparms') can only include: ",paste(parmlist[-length(parmlist)],collapse=", "),", or ",parmlist[length(parmlist)]))
  
  commonparms <- parms
  
  getlogitp <- derivedlogitfun(parms,"p")
  getlogitc <- derivedlogitfun(parms,"c")  
  
  if(any(parms=="p")){
    namesp <- paste0("p[",1:noccas,"]")   
    commonparms <- commonparms[-which(parms=="p")]
    parms <- parms[-which(parms=="p")]
    parms <- c(parms,namesp)
  } else {
    namesp <- NULL
  }
  if(any(parms=="c")){
    namesc <- paste0("c[",2:noccas,"]")   
    commonparms <- commonparms[-which(parms=="c")]
    parms <- parms[-which(parms=="c")]
    parms <- c(parms,namesc)
  } else {
    namesc <- NULL
  }
  list(commonparms=commonparms,parms=parms,namesp=namesp,namesc=namesc,getlogitp=getlogitp,getlogitc=getlogitc)
}

checkmmClosedinput<-function(mmslist,modlist,nmod,nchains,iter,miter,mburnin,mthin,modprior,M1){
  if(!all(match(unlist(unique(lapply(modlist,names))),c("mcmc","mod.p","mod.delta","DM","initial.values","priorparms","mms"),nomatch=0))) stop("each object in 'modlist' must be a list returned by multimarkClosed() or markClosed()")
  if(!all(lapply(modlist,function(x) is.mcmc.list(x$mcmc))==TRUE)) stop("mcmc output for each model must be an object of type 'mcmc.list'")
  if(nmod<2) stop("'modlist' must contain at least two models")
  if(length(mmslist)!=1) stop("'multimarksetup' (mms) object must be identical for each model")
  if(length(nchains)!=1) stop("all models must have same number of chains")
  if(length(iter)!=1) stop("all chains must have same number of iterations")
  if(miter<=mburnin) stop("'mburnin' must be less than ",miter)
  if(mthin>max(1,floor((miter-mburnin+1)/2)) | mthin<1) stop("'mthin' must be >0 and <=",max(1,floor((miter-mburnin+1)/2)))
  if(length(modprior)!=nmod | base::sum(modprior)!=1) stop(paste("'modprior' must be a vector of length ",nmod," that sums to 1"))
  if(length(M1)!=nchains) stop("'M1' must be an integer vector of length ",nchains)
  if(!all(match(M1,1:nmod,nomatch=0))) stop("'M1' must be an integer vector of length ",nchains," with values ranging from 1 to ",nmod)
  mms<-mmslist[[1]]
  if(class(mms)!="multimarksetup") stop("'mms' for each model must be an object of class 'multimarksetup'")
  return(mms)
}

drawmissingClosed<-function(M.cur,missing,pbetapropsd,sigppropshape,sigppropscale){
  missingpbeta <- rnorm(length(missing$missingpbetaparms[[M.cur]]),sd=pbetapropsd)
  names(missingpbeta) <- missing$missingpbetaparms[[M.cur]]
  missingdelta <- numeric()
  if(length(missing$missingdeltaparms[[M.cur]])){
    if(length(missing$missingdeltaparms[[M.cur]])==1){
      missingdelta <- rbeta(1,1,1)/2
    } else {
      missingdelta <- rdirichlet(1,c(1,1,1))[1:2]
    }
  }
  names(missingdelta) <- missing$missingdeltaparms[[M.cur]]
  missingsigp  <- rinvgamma(length(missing$missingsigpparms[[M.cur]]),shape=sigppropshape,scale=sigppropscale)
  names(missingsigp) <- missing$missingsigpparms[[M.cur]]
  missingzp <- rnorm(length(missing$missingzpparms[[M.cur]]),sd=missing$zppropsd+sqrt(missingsigp)*missing$usesigp)
  names(missingzp) <- missing$missingzpparms[[M.cur]]
  missing <- c(missingpbeta,missingdelta,missingzp,missingsigp)
  missing
}

getbrobprobClosed<-function(imod,modprior,posterior,cur.parms,missing,pbetapropsd,sigppropshape,sigppropscale){
  deltadens <- 0
  if(length(missing$missingdeltaparms[[imod]])){
    if(length(missing$missingdeltaparms[[imod]])==1){
      deltadens <- dbeta(2*cur.parms[missing$missingdeltaparms[[imod]]],1,1,log=TRUE)
    } else {
      delta <- c(cur.parms[missing$missingdeltaparms[[imod]]],1.-sum(cur.parms[missing$missingdeltaparms[[imod]]]))
      deltadens <- ddirichlet(delta,c(1,1,1))
    }
  }
  brob(log(modprior[imod]) 
       + posterior 
       + base::sum(dnorm(cur.parms[missing$missingpbetaparms[[imod]]],sd=pbetapropsd,log=TRUE))
       + base::sum(dnorm(cur.parms[missing$missingzpparms[[imod]]],sd=missing$zppropsd+sqrt(cur.parms[missing$missingsigpparms[[imod]]])*missing$usesigp,log=TRUE))
       + base::sum(dinvgamma(cur.parms[missing$missingsigpparms[[imod]]],shape=sigppropshape,scale=sigppropscale))
       + deltadens)
}

getcurClosedparmslist<-function(cur.parms,DM,M,noccas,data_type,alpha){
  
  parmslist=vector('list',1)
  parmslist[[1]]$H<-cur.parms[paste0("H[",1:M,"]")]
  parmslist[[1]]$N <- cur.parms["N"]
  parmslist[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
  parmslist[[1]]$zp <- cur.parms[paste0("zp[",1:M,"]")]
  parmslist[[1]]$sigma2_zp <- cur.parms["sigma2_zp"]
  
  parmslist[[1]]$psi <- cur.parms["psi"]
  parmslist[[1]]$delta_1 <- cur.parms["delta_1"]
  parmslist[[1]]$delta_2 <- cur.parms["delta_2"]
  parmslist[[1]]$delta <- cur.parms["delta"]

  if(data_type=="sometimes"){
    parmslist[[1]]$alpha <- cur.parms["alpha"]
  } else {
    parmslist[[1]]$alpha <- alpha   
  }
  parmslist
}

missingparmnamesClosed<-function(params,M,noccas,zppropsd){
  
  multiparms <- unique(unlist(params))
  
  commonparms <- Reduce(intersect, params)
  commonparms <- commonparms[-match("logPosterior",commonparms)]
  
  missingparms <- lapply(params,get_missingparms,multiparms=multiparms)
  
  missingpbetaparms <- extractmissingparms(missingparms,"pbeta")
  
  missingdeltaparms <- extractmissingparms(missingparms,"delta")
  
  missingsigpparms <- lapply(missingparms,function(x) unlist(x,use.names=FALSE)[which(x=="sigma2_zp")])
  
  missingzpparms <- extractmissingparms(missingparms,"zp[")
  
  if(is.null(zppropsd)){
    zppropsd <- 0
    usesigp <- 1
  } else {
    usesigp <-0
  }
  list(commonparms=commonparms,missingparms=missingparms,missingpbetaparms=missingpbetaparms,missingdeltaparms=missingdeltaparms,missingsigpparms=missingsigpparms,missingzpparms=missingzpparms,zppropsd=zppropsd,usesigp=usesigp) 
}

rjmcmcClosed <- function(ichain,mms,M,noccas,data_type,alpha,C,All.hists,modlist,DMlist,deltalist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,sigppropshape,sigppropscale,pmodnames,deltamodnames,gq,printlog){
  
  multimodel <- matrix(0,nrow=(max(1,floor(miter/mthin)))-(floor(mburnin/mthin)),ncol=length(monitorparms$parms)+1,dimnames=list(NULL,c(monitorparms$parms,"M")))
  
  nmod <- length(modlist)
  mod.prob.brob <- as.brob(numeric(nmod))
  
  commonparms <- monitorparms$commonparms
  
  if(any(deltalist==~NULL)){
    H<-get_H(mms,mms@naivex)
    names(H)<-paste0("H[",1:M,"]")
  } else {
    H<-NULL
  }
  
  M.cur<- M1
  
  modmissingparms <- drawmissingClosed(M.cur,missing,pbetapropsd,sigppropshape,sigppropscale)
  cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
  
  DM <- DMlist[[M.cur]]
  DM$mod.delta <- deltalist[[M.cur]]
  DM$mod.p.h <- mod.p.h[[M.cur]]
  
  cur.parms.list <- getcurClosedparmslist(cur.parms,DM,M,noccas,data_type,alpha)  
  
  for(iiter in 1:miter){
    
    mod.prob.brob[M.cur] <- getbrobprobClosed(M.cur,modprior,cur.parms["logPosterior"],cur.parms,missing,pbetapropsd,sigppropshape,sigppropscale)
    
    for(imod in (1:nmod)[-M.cur]){ 
      
      DM <- DMlist[[imod]]
      DM$mod.delta <- deltalist[[imod]]
      DM$mod.p.h <- mod.p.h[[imod]]
      
      cur.parms.list[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
      
      loglike <- loglikeClosed(cur.parms.list[[1]],DM,noccas,C,All.hists,gq[[imod]])
      
      posterior <- loglike + priorsClosed(cur.parms.list[[1]],DM,priorlist[[imod]],data_type)
      
      mod.prob.brob[imod] <- getbrobprobClosed(imod,modprior,posterior,cur.parms,missing,pbetapropsd,sigppropshape,sigppropscale)
    }
    
    if(any(is.na(as.numeric(mod.prob.brob)))){
      warning(paste0("'NA' posterior for model '","p(",pmodnames[is.na(as.numeric(mod.prob.brob))],")delta(",deltamodnames[is.na(as.numeric(mod.prob.brob))],")' at iteration ",iiter,"; model move rejected."))
      flush.console()
    } else {       
      mod.prob <- as.numeric(mod.prob.brob/Brobdingnag::sum(mod.prob.brob))
      M.cur <- (1:nmod)[rmultinom(1, 1, mod.prob)==1]
    }
    
    modmissingparms <- drawmissingClosed(M.cur,missing,pbetapropsd,sigppropshape,sigppropscale)
    cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
    
    DM <- DMlist[[M.cur]]
    DM$mod.delta <- deltalist[[M.cur]]
    DM$mod.p.h <- mod.p.h[[M.cur]]
    
    cur.parms.list <- getcurClosedparmslist(cur.parms,DM,M,noccas,data_type,alpha)  
    
    if(iiter>mburnin & !iiter%%mthin){
      multimodel[iiter/mthin-floor(mburnin/mthin),"M"] <- M.cur
      multimodel[iiter/mthin-floor(mburnin/mthin),commonparms] <- cur.parms[commonparms]
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesp] <- monitorparms$getlogitp(DM$mod.p.h,DM$p,cur.parms.list[[1]]$pbeta,cur.parms.list[[1]]$sigma2_zp)
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesc] <- monitorparms$getlogitc(DM$mod.p.h,DM$c,cur.parms.list[[1]]$pbeta,cur.parms.list[[1]]$sigma2_zp)[-1]
    }
    
    if(!(iiter%%(miter/ min(miter,100)))) {
      if(printlog){
        cat("Chain ",ichain," is ",100*(iiter/miter),"% complete \n",sep="")        
      } else{
        cat("\rChain ",ichain," is ",100*(iiter/miter),"% complete",sep="")
      }
    }
  }
  return(multimodel)
}

#' Multimodel inference for 'multimark' closed population abundance models
#'
#' This function performs Bayesian multimodel inference for a set of 'multimark' closed population abundance models using the reversible jump Markov chain Monte Carlo (RJMCMC) algorithm proposed by Barker & Link (2013).
#'
#'
#' @param modlist A list of individual model output lists returned by \code{\link{multimarkClosed}}. The models must have the same number of chains and MCMC iterations.
#' @param modprior Vector of length \code{length(modlist)} containing prior model probabilities. Default is \code{modprior = rep(1/length(modlist), length(modlist))}.
#' @param monparms Parameters to monitor. Only parameters common to all models can be monitored (e.g., "\code{pbeta[(Intercept)]}", "\code{N}", "\code{psi}"), but derived capture ("\code{p}") and recapture ("\code{c}") probabilities can also be monitored. Default is \code{monparms = "N"}.
#' @param miter The number of RJMCMC iterations per chain. If \code{NULL}, then the number of MCMC iterations for each individual model chain is used.
#' @param mburnin Number of burn-in iterations (\code{0 <= mburnin < miter}).
#' @param mthin Thinning interval for monitored parameters.
#' @param M1 Integer vector indicating the initial model for each chain, where \code{M1_j=i} initializes the RJMCMC algorithm for chain j in the model corresponding to \code{modlist[[i]]} for i=1,...,  \code{length(modlist)}. If \code{NULL}, the algorithm for all chains is initialized in the most general model. Default is \code{M1=NULL}.
#' @param pbetapropsd Scaler specifying the standard deviation of the Normal(0, pbetapropsd) proposal distribution for "\code{pbeta}"  parameters. Default is \code{pbetapropsd=1}. See Barker & Link (2013) for more details.
#' @param zppropsd Scaler specifying the standard deviation of the Normal(0, zppropsd) proposal distribution for "\code{zp}"  parameters. Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. If \code{NULL}, zppropsd = sqrt(sigma2_zp) is used. Default is \code{zppropsd=NULL}. See Barker & Link (2013) for more details.  
#' @param sigppropshape Scaler specifying the shape parameter of the invGamma(shape = sigppropshape, scale = sigppropscale) proposal distribution for "\code{sigma_zp=sqrt(sigma2_zp)}". Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. Default is \code{sigppropshape=6}. See Barker & Link (2013) for more details.
#' @param sigppropscale Scaler specifying the scale parameter of the invGamma(shape = sigppropshape, scale = sigppropscale) proposal distribution for "\code{sigma_zp=sqrt(sigma2_zp)}". Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. Default is \code{sigppropscale=4}. See Barker & Link (2013) for more details.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @details Note that setting \code{parms="all"} is required when fitting individual \code{\link{multimarkClosed}} models to be included in \code{modlist}.
#' @return A list containing the following:
#' \item{rjmcmc}{Reversible jump Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}. Includes RJMCMC output for monitored parameters and the current model at each iteration ("\code{M}").}
#' \item{pos.prob}{A list of calculated posterior model probabilities for each chain, including the overall posterior model probabilities across all chains.}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkClosed}}, \code{\link{processdata}}
#' @references
#' Barker, R. J. and Link. W. A. 2013. Bayesian multimodel inference by RJMCMC: a Gibbs sampling approach. The American Statistician 67: 150-156.
#' @examples
#' \dontshow{
#' setup<-processdata(bobcat)
#' test.dot<-multimarkClosed(mms=setup,parms="all",iter=10,burnin=0,bin=5)
#' test<-multimodelClosed(modlist=list(mod1=test.dot,mod2=test.dot))
#' }
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Generate object of class "multimarksetup"
#' setup <- processdata(bobcat)
#'  
#' #Run single chain using the default model for bobcat data. Note parms="all".
#' bobcat.dot <- multimarkClosed(mms=setup,parms="all")
#' 
#' #Run single chain for bobcat data with time effects. Note parms="all".
#' bobcat.time <- multimarkClosed(mms=setup,mod.p=~time,parms="all")
#' 
#' #Perform RJMCMC using defaults
#' modlist <- list(mod1=bobcat.dot,mod2=bobcat.time)
#' bobcat.M <- multimodelClosed(modlist=modlist,monparms=c("N","p"))
#' 
#' #Posterior model probabilities
#' bobcat.M$pos.prob
#'  
#' #multimodel posterior summary for abundance
#' summary(bobcat.M$rjmcmc[,"N"])}
multimodelClosed<-function(modlist,modprior=rep(1/length(modlist),length(modlist)),monparms="N",miter=NULL,mburnin=0,mthin=1,M1=NULL,pbetapropsd=1,zppropsd=NULL,sigppropshape=6,sigppropscale=4,printlog=FALSE){
  
  nmod <- length(modlist)
  iter <- unlist(unique(lapply(modlist,function(x) unique(lapply(x$mcmc,nrow)))))
  nchains <- unlist(unique(lapply(modlist,function(x) length(x$mcmc))))
  mmslist <- unlist(unique(lapply(modlist, function(x) {x$mms@covs<-data.frame();x$mms})))
  
  params <- lapply(modlist,function(x) varnames(x$mcmc))
  
  if(is.null(M1)) M1 <- rep(which.max(lapply(params,length))[1],nchains)
  
  if(is.null(miter)) miter <- iter
  
  mms<-checkmmClosedinput(mmslist,modlist,nmod,nchains,iter,miter,mburnin,mthin,modprior,M1)
  
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas)
  C<-mms@C
  gq <- lapply(modlist,function(x) gauss.quad(x$priorparms$npoints,kind="hermite"))
  
  checkparmsClosed(mms,modlist,params,parmlist=c("pbeta[(Intercept)]","N","logPosterior"),M)
  
  pmodnames <- unlist(lapply(modlist,function(x) x$mod.p)) 
  deltamodnames <- unlist(lapply(modlist,function(x) x$mod.delta)) 
  
  message("\nPerforming closed population Bayesian multimodel inference by RJMCMC \n")
  if(all(deltamodnames!=~NULL)) {
    message(paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,")\n"))
  } else if(all(deltamodnames==~NULL)){
    message(paste0("mod",1:nmod,": ","p(",pmodnames,")\n"))
  }
  
  missing <- missingparmnamesClosed(params,M,noccas,zppropsd) 
  
  monitorparms <- monitorparmsClosed(monparms,c(missing$commonparms,"p","c"),noccas)
  
  DMlist <- lapply(modlist,function(x) x$DM)
  deltalist <- lapply(modlist,function(x) x$mod.delta)
  priorlist <- lapply(modlist,function(x) x$priorparms) 
  mod.p.h <- unlist(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.p))$term.labels)))
  
  data_type <- mms@data.type
  if(data_type=="never"){
    alpha <- 0
  } else if(data_type=="always"){
    alpha <- 1
  } else {
    alpha <- numeric(0)
  }
 

  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0("multimodelClosed_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("rjmcmcClosed"),envir=environment())
    multimodel <- parLapply(cl,1:nchains, function(ichain) rjmcmcClosed(ichain,mms,M,noccas,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[ichain]]),DMlist,deltalist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1[ichain],monitorparms,missing,pbetapropsd,sigppropshape,sigppropscale,pmodnames,deltamodnames,gq,printlog))
    stopCluster(cl)
    gc()
  } else {
    multimodel <- vector('list',nchains)
    multimodel[[nchains]] <- rjmcmcClosed(nchains,mms,M,noccas,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[nchains]]),DMlist,deltalist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,sigppropshape,sigppropscale,pmodnames,deltamodnames,gq,printlog)
    gc()
  }
  
  if(mburnin<mthin){
    temp=seq(mthin,max(1,miter),mthin)
  } else {
    temp=seq(mthin*(floor(mburnin/mthin)+1),miter,mthin)
  }
  
  pos.prob <- vector('list',nchains)
  for(ichain in 1:nchains){
    pos.prob[[ichain]] <-hist(multimodel[[ichain]][,"M"],plot=F,breaks=0:nmod)$density
    if(all(deltamodnames!=~NULL)){
      names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,")") 
    } else {
      names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")")
    }
    multimodel[[ichain]] <- mcmc(multimodel[[ichain]])
    attributes(multimodel[[ichain]])$mcpar <- c(head(temp,n=1),tail(temp,n=1),mthin)
  }  
  
  multimodel <- as.mcmc.list(multimodel)
  names(pos.prob) <- paste0("chain",1:nchains)
  pos.prob[["overall"]]<- hist(unlist(multimodel[, "M"]),plot = F, breaks = 0:nmod)$density
  if(all(deltamodnames!=~NULL)){
    names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,")") 
  } else {
    names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")")
  }
  list(rjmcmc=multimodel,pos.prob=pos.prob) 
}