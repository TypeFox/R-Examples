#' Simulate open population capture-mark-recapture data arising from multiple non-invasive marks
#' 
#' This function generates encounter histories from simulated open population capture-mark-recapture data consisting of multiple non-invasive marks. 
#' 
#' 
#' @param N Number of individuals.
#' @param noccas Number of sampling occasions. \code{floor(N/noccas)} individuals are first encountered on each occasion.
#' @param pbeta Logit- or probit-scale intercept term(s) for capture probability (p). Must be a scaler or vector of length \code{noccas}.
#' @param sigma2_zp Logit- or probit-scale individual heterogeneity variance term for capture probability (p).
#' @param phibeta Logit- or probit-scale intercept term(s) for survival probability (\eqn{\phi}). Must be a scaler or vector of length \code{noccas}.
#' @param sigma2_zphi Logit- or probit-scale individual heterogeneity variance term for survival probability (\eqn{\phi}).
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
#' @param link Link function for detection probability. Must be "\code{logit}" or "\code{probit}". Note that \code{\link{multimarkCJS}} is currently implemented for the probit link only.
#'
#' @return A list containing the following:
#' \item{Enc.Mat}{A matrix containing the observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions.}
#' \item{trueEnc.Mat}{A matrix containing the true (latent) encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions.}
#' @author Brett T. McClintock 
#' @seealso \code{\link{processdata}}, \code{\link{multimarkCJS}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' #simulate data for data.type="sometimes" using defaults
#' data<-simdataCJS(data.type="sometimes")
simdataCJS <- function(N=100,noccas=5,pbeta=-0.25,sigma2_zp=0,phibeta=1.6,sigma2_zphi=0,delta_1=0.4,delta_2=0.4,alpha=0.5,data.type="never",link="probit"){
  
  if(length(pbeta)==1){
    pbeta=rep(pbeta,noccas)
  } else if(length(pbeta)!=noccas){
    stop(paste0("'pbeta' must be a scaler or vector of length ",noccas))
  }
  if(length(phibeta)==1){
    phibeta=rep(phibeta,noccas-1)
  } else if(length(phibeta)!=noccas-1){
    stop(paste0("'phibeta' must be a scaler or vector of length ",noccas-1))
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
  
  n <- floor(N/noccas)
  tEnc.Mat <-matrix(0,nrow=N,ncol=noccas)
  first <- sort(rep(1:noccas,n))
  zp <- rnorm(N,0,sqrt(sigma2_zp))
  zphi <- rnorm(N,0,sqrt(sigma2_zphi))
  tmp.phibeta <- c(phibeta,0)
  for(i in 1:(n*(noccas-1))){
    ind <- tEnc.Mat[i,first[i]] <- 1
    for(j in (first[i]+1):noccas){
      if(link=="probit"){
        p<-pnorm(pbeta[j]+zp[i])
        phi<-pnorm(tmp.phibeta[j]+zphi[i])
      } else if(link=="logit"){
        p<-expit(pbeta[j]+zp[i])      
        phi<-expit(tmp.phibeta[j]+zphi[i])  
      } else {stop("link function must be 'probit' or 'logit'")}
      tEnc.Mat[i,j]<-rbinom(1,1,p*ind)
      if(runif(1)>phi){
        ind<-0
      }
    }
  }
  for(i in ((noccas-1)*n + 1):N){
    tEnc.Mat[i,first[i]]<-1
  }
  Rand.Mat<-matrix(runif(N*noccas,0,1),N,noccas)
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat<delta_2)] <- 2      # type 2 encounters
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat>(1-delta_B))] <- 4  # type 1 and type 2 encounters
  tEnc.Mat[which(tEnc.Mat==4)] <- tEnc.Mat[which(tEnc.Mat==4)]-(runif(base::sum(tEnc.Mat==4))<(1-alpha))   # unobserved type 1 and type 2 encounters
  Enc.Mat <- get_Enc(tEnc.Mat,data.type)
  return(list(Enc.Mat=Enc.Mat,trueEnc.Mat=tEnc.Mat))
}

get_q<-function(mms,DM,H,pbeta,zp,phibeta,zphi){
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  
  p <- array(0,dim=c(noccas-1,M,noccas-1))
  phi <- array(0,dim=c(noccas-1,M,noccas-1))  
  
  XBp=DM$p%*%pbeta
  XBphi=DM$phi%*%phibeta   
  
  for(cohort in 1:(noccas-1)){
    ind <- c(0,cumsum(seq(noccas-1,2)))[cohort]+1:(noccas-cohort)
    p[cohort,,cohort:(noccas-1)] <- pmin(pmax(pnorm(matrix(rep(XBp[ind],each=M)+zp,nrow=M,ncol=noccas-cohort)),tol),1.-tol) 
    phi[cohort,,cohort:(noccas-1)] <- pmin(pmax(pnorm(matrix(rep(XBphi[ind],each=M)+zphi,nrow=M,ncol=noccas-cohort)),tol),1.-tol)  
  }
  
  q <- matrix(0,M,noccas)
  
  for(i in 1:M){
    firstcap <- mms@C[H[i]]
    lastcap <- mms@L[H[i]]
    if(H[i]>1){
      q[i,firstcap:lastcap] <- 1
      if(lastcap<noccas){
        for(t in (lastcap+1):noccas){
          num <- phi[firstcap,i,t-1] * q[i,t-1] * (1.-p[firstcap,i,t-1]);
          if(t<noccas) num <- num * (1.-phi[firstcap,i,t])
          probq <- num/(num + (1.-phi[firstcap,i,t-1]*q[i,t-1]));
          if(t<noccas){
            if(q[i,t+1]) probq <- 1
          }
          q[i,t] <- rbinom(1,1.0,probq);
        }
      }
    }
  }
  q
}

#sampleZ<-function(mms,DM,Hs,pbeta,zp,phibeta,zphi){
#  noccas<-ncol(mms@Enc.Mat)
#  M<-nrow(mms@Enc.Mat)
#  
#  Allhists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas)
#  
#  p <- array(0,dim=c(noccas-1,M,noccas-1))
#  phi <- array(0,dim=c(noccas-1,M,noccas-1))  
#  
#  XBp=DM$p%*%pbeta
#  XBphi=DM$phi%*%phibeta   
#  
#  for(cohort in 1:(noccas-1)){
#    ind <- c(0,cumsum(seq(noccas-1,2)))[cohort]+1:(noccas-cohort)
#    p[cohort,,cohort:(noccas-1)] <- pmin(pmax(pnorm(matrix(rep(XBp[ind],each=M)+zp,nrow=M,ncol=noccas-cohort)),tol),1.-tol)  
#    phi[cohort,,cohort:(noccas-1)] <- pmin(pmax(pnorm(matrix(rep(XBphi[ind],each=M)+zphi,nrow=M,ncol=noccas-cohort)),tol),1.-tol)  
#  }
#  
#  q <- matrix(0,M,noccas)
#  
#  firstcap <- mms@C[Hs]
#  
#  pi<-array(0,dim=c(M,noccas,2));
#  
#  prob=numeric(2);
#
#  for(i in 1:M){
#    if(Hs[i]>1){
#      q[i,firstcap[i]] = 1;
#      if(firstcap[i]<noccas){
#        indhist= Allhists[Hs[i],firstcap[i]+1];
#        if(indhist) {
#          pi[i,firstcap[i]+1,1] = 0.;
#          pi[i,firstcap[i]+1,2] = 1.;
#        } else {
#          pi[i,firstcap[i]+1,1] = 1. - phi[firstcap[i],i,firstcap[i]];
#          pi[i,firstcap[i]+1,2] = phi[firstcap[i],i,firstcap[i]] * (1.-p[firstcap[i],i,firstcap[i]]);
#        }
#        if(firstcap[i]+1<noccas){
#          for(t in (firstcap[i]+2):noccas){
#            indhist= Allhists[Hs[i],t];
#            if(indhist) {
#              pi[i,t,1] = 0.;
#              pi[i,t,2] = 1.;
#            } else {
#              pi[i,t,1] = pi[i,t-1,1] + pi[i,t-1,2]*(1.-phi[firstcap[i],i,t-2]);
#              pi[i,t,2] = pi[i,t-1,2] * phi[firstcap[i],i,t-2] *(1.-p[firstcap[i],i,t-1]); 
#            } 
#          }
#        }
#        q[i,noccas] = sample(0:1,1,prob=pi[i,noccas,]);
#        if(firstcap[i]<(noccas-1)){
#          for(t in (noccas-1):(firstcap[i]+1)){
#            prob[1] = pi[i,t,1] * (1.-q[i,t+1]);
#            prob[2] = pi[i,t,2] * (q[i,t+1]*phi[firstcap[i],i,t]*(1.-p[firstcap[i],i,t])  + (1.-q[i,t+1])*(1.-phi[firstcap[i],i,t]));
#            q[i,t] = sample(0:1,1,prob=prob);
#          }
#        }     
#      }
#    }
#  }
#  q
#}  

getDMformula<-function(mod){
  mod.terms<-attributes(terms(mod))$term.labels
  if(any(mod.terms=="c")) stop("'c' is not allowed in model formulae for 'CJS' models")
  if(any(mod.terms=="time") | any(mod.terms=="age")){
    mod<-formula(paste("~-1",paste(mod.terms,collapse="+"),sep="+"))
    mod.terms<-attributes(terms(mod))$term.labels
  }
  if(any(mod.terms=="h")){
    if(length(mod.terms)>1){
      mod<-formula(paste("~-1",paste(mod.terms[-which(mod.terms=="h")],collapse="+"),sep="+"))
    } else {
      mod<-formula(~1)  
    }
    mod.h<-TRUE 
  } else {
    mod.h<-FALSE
  }
  mod.terms<-attributes(terms(mod))$term.labels
  list(mod=mod,mod.terms=mod.terms,mod.h=mod.h)
}

formatDM<-function(mod,temp,noccas,parm){
  DM <- model.matrix(mod$mod,temp)
  modterms <- mod$mod.terms
  if(length(modterms)){
    for(i in 1:length(modterms)){
      if(modterms[i]=="time"){
        if(length(unique(temp[,modterms[i]]))!=sum(attr( DM,"assign")==i)){
          DM<-DM[,-which(attr( DM,"assign")==i)]
          DM<-cbind(DM,model.matrix(~-1+time,temp))
        }
      }
      if(modterms[i]=="age"){
        if(length(unique(temp[,modterms[i]]))!=sum(attr( DM,"assign")==i)){
          DM<-DM[,-which(attr( DM,"assign")==i)]
          DM<-cbind(DM,model.matrix(~-1+age,temp))
        }
      }
    }
  }
  if(nrow(DM)!=sum((noccas-1):1)) stop(paste("model design matrix for",parm,"p must have",sum((noccas-1):1),"rows"))
  if(!any(colnames(DM)=="(Intercept)")){
    DM<-cbind(rep(1,sum((noccas-1):1)),DM)
    colnames(DM)[1]<-"(Intercept)"
  }
  DM
}

get_DMCJS<-function(mod.p,mod.phi,mod.delta,Enc.Mat,covs,type="CJS",...){
  Enc.Mat[which(Enc.Mat>0)] <- 1
  ch<-as.character(as.matrix( apply(Enc.Mat, 1, paste, collapse = ""), ncol=1 ))
  if(length(which(ch==paste(rep(0,ncol(Enc.Mat)),collapse="")))){
    ch<-ch[-which(ch==paste(rep(0,ncol(Enc.Mat)),collapse=""))] 
  }
  CH<-process.data(data.frame(ch,options(stringsAsFactors = FALSE)),model=type)
  temp<-make.design.data(CH,...)
  if((CH$nocc-1) %% 2 == 0){
    temp$Phi$Time <- temp$Phi$Time+1-((CH$nocc-1)/2)-.5
    temp$p$Time <- temp$p$Time+1-((CH$nocc-1)/2)-.5
  } else {
    temp$Phi$Time <- temp$Phi$Time+1-ceiling((CH$nocc-1)/2)
    temp$p$Time <- temp$p$Time+1-ceiling((CH$nocc-1)/2)
  }
  
  pmod <- getDMformula(mod.p)
  phimod <- getDMformula(mod.phi)
  
  if(length(covs)){
    covs <- covs[-1,,drop=FALSE]
    temp$p<-cbind(temp$p,covs[temp$p$time,,drop=FALSE])
    rownames(temp$p)<-NULL
    temp$Phi<-cbind(temp$Phi,covs[temp$Phi$time,,drop=FALSE])
    rownames(temp$Phi)<-NULL
  }
  
  DMp <-formatDM(pmod,temp$p,CH$nocc,"p")
  DMphi <-formatDM(phimod,temp$Phi,CH$nocc,"phi")  
  
  pphiseq<-matrix(0,(CH$nocc-1),(CH$nocc-1))
  for(t in 1:(CH$nocc-1)){
    pphiseq[t:(CH$nocc-1),t]<-t:(CH$nocc-1)
  }
  rownames(DMp) <- paste0("p[",rep(1:(CH$nocc-1),times=(CH$nocc-1):1),",",pphiseq[which(pphiseq>0)]+1,"]")
  rownames(DMphi) <- paste0("phi[",rep(1:(CH$nocc-1),times=(CH$nocc-1):1),",",pphiseq[which(pphiseq>0)],"]")
  
  deltattr <- attributes(terms(mod.delta))$term.labels
  if(length(deltattr)){
    if(deltattr!="type") stop("'mod.delta' must be '~1' or '~type'")
  }
  return(list(p=DMp,mod.p.h=pmod$mod.h,phi=DMphi,mod.phi.h=phimod$mod.h,mod.delta=mod.delta))
}

mcmcCJS<-function(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sdp,Prop.sdphi,maxnumbasis,pbeta0,pprec0,phibeta0,phiprec0,l0p,d0p,l0phi,d0phi,a0delta,a0alpha,b0alpha,a0psi,b0psi,link,printlog){
  
  M<-nrow(mms@Enc.Mat)
  noccas<-ncol(mms@Enc.Mat)
  DMp<-DM$p
  DMphi<-DM$phi
  mod.p.h<-DM$mod.p.h
  pdim<-ncol(DMp)
  mod.phi.h<-DM$mod.phi.h
  phidim<-ncol(DMphi)
  
  #declare and initialize parameters
  pbeta<-rep(NA,max(1,floor(iter/thin))*(pdim))
  zp<-rep(NA,ifelse(any(params=="zp"),max(1,floor(iter/thin))*M,M))
  phibeta<-rep(NA,max(1,floor(iter/thin))*(phidim))
  zphi<-rep(NA,ifelse(any(params=="zphi"),max(1,floor(iter/thin))*M,M))
  H<-rep(NA,ifelse(any(params=="H"),max(1,floor(iter/thin))*M,M))
  q<-rep(NA,ifelse(any(params=="q"),max(1,floor(iter/thin))*M*noccas,M*noccas))
  sigma2_zp<-rep(NA,max(1,floor(iter/thin)))
  sigma2_zphi<-rep(NA,max(1,floor(iter/thin)))
  alpha<-rep(NA,max(1,floor(iter/thin)))
  delta_1<-rep(NA,max(1,floor(iter/thin)))
  delta_2<-rep(NA,max(1,floor(iter/thin)))
  psi<-rep(NA,max(1,floor(iter/thin)))
  loglike<-rep(NA,max(1,floor(iter/thin)))
  
  pbeta[1:pdim] <- inits[[ichain]]$pbeta
  phibeta[1:phidim] <- inits[[ichain]]$phibeta
  zp[1:M] <- inits[[ichain]]$zp
  zphi[1:M] <- inits[[ichain]]$zphi
  H[1:M] <- inits[[ichain]]$H-1
  q[1:(noccas*M)] <- c(t(inits[[ichain]]$q))
  sigma2_zp[1] <- inits[[ichain]]$sigma2_zp
  sigma2_zphi[1] <- inits[[ichain]]$sigma2_zphi
  alpha[1] <- inits[[ichain]]$alpha
  delta_1[1] <- inits[[ichain]]$delta_1
  delta_2[1] <- inits[[ichain]]$delta_2
  psi[1] <- inits[[ichain]]$psi
  
  if(link=="probit"){
    
    posterior <- .C(ProbitCJSC,as.integer(ichain),as.numeric(pbeta0),as.numeric(c(pprec0)),as.numeric(pbeta),as.numeric(phibeta0),as.numeric(c(phiprec0)),as.numeric(phibeta),as.numeric(zp),as.numeric(sigma2_zp),as.numeric(zphi),as.numeric(sigma2_zphi),as.numeric(delta_1),as.numeric(delta_2),as.numeric(alpha),as.integer(inits[[ichain]]$x),as.numeric(psi), as.integer(H), as.integer(q),
                    as.integer(noccas), as.integer(M), as.numeric(a0delta), as.numeric(a0alpha), as.numeric(b0alpha),as.numeric(l0p),as.numeric(d0p),as.numeric(l0phi),as.numeric(d0phi),as.numeric(a0psi),as.numeric(b0psi),
                    as.numeric(loglike),
                    as.integer(length(mms@vAll.hists)/noccas),as.integer(mms@vAll.hists),as.integer(mms@C),as.integer(mms@L),as.integer(mms@indBasis-1), as.integer(mms@ncolbasis), as.integer(mms@knownx),as.numeric(as.vector(t(DMp))),as.numeric(as.vector(t(DMphi))),as.integer(dim(DMp)),as.integer(dim(DMphi)),
                    as.integer(iter), as.integer(thin),as.integer(maxnumbasis),
                    as.integer(mod.p.h),as.integer(mod.phi.h),as.integer(mms@data.type=="sometimes"),as.integer(any(params=="zp")),as.integer(any(params=="zphi")),as.integer(any(params=="q")),as.integer(any(params=="H")),as.integer(DM$mod.delta != ~NULL),as.integer(DM$mod.delta==formula(~type)),as.integer(printlog),NAOK = TRUE)
  } else {
    stop("only 'probit' link is currently implemented for CJS models")
  }
  names(posterior) <- c("ichain","pbeta0","pprec0","pbeta","phibeta0","phiprec0","phibeta", "zp", "sigma2_zp", "zphi", "sigma2_zphi", "delta_1","delta_2","alpha", "x","psi","H","q",
                        "noccas", "M","a0delta", "a0alpha", "b0alpha","l0p","d0p","l0phi","d0phi","a0psi","b0psi",
                        "loglike",
                        "nHists","vAll.hists","C","L", "indBasis", "ncolBasis","knownx","DMp","DMphi","pdim","phidim",
                        "iter", "thin", "maxnumbasis",
                        "mod.p.h","mod.phi.h","sometimes?","zp?","zphi?","q?","H?","updatedelta?","type","printlog?")
  
  g <- posterior$iter
  x <- posterior$x
  if(any(params=="zp")){
    zp<-matrix(posterior$zp[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)],ncol=M,byrow=T)
    zpout<-NULL
  } else {
    zp <- NULL
    zpout <- posterior$zp
  }
  temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),
              matrix(posterior$phibeta[(floor(burnin/thin)*phidim+1):(max(1,floor(iter/thin))*phidim)],ncol=phidim,byrow=T),
              posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],
              zp,posterior$sigma2_zp[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
  if(any(params=="zphi")){
    zphi <- matrix(posterior$zphi[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)],ncol=M,byrow=T)
    zphiout <- NULL
  } else {
    zphi <- NULL
    zphiout <- posterior$zphi
  }
  temp <- cbind(temp,zphi,posterior$sigma2_zphi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
  if(any(params=="q")){
    q<-matrix(posterior$q[(floor(burnin/thin)*M*noccas+1):(max(1,floor(iter/thin))*M*noccas)],ncol=M*noccas,byrow=T)
    qout <- NULL
  } else {
    q <- NULL
    qout <- posterior$q
  }
  temp<-cbind(temp,q) 
  if(any(params=="H")){
    H <- matrix(posterior$H[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)]+1,ncol=M,byrow=T)
    Hout <- NULL
  } else {
    H <- NULL
    Hout <- posterior$H+1
  }
  posterior<-cbind(temp,H,posterior$loglike[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
  return(list(posterior=posterior,x=x,H=Hout,q=qout,zp=zpout,zphi=zphiout,g=g))
}

loglikeCJS<-function(parms,DM,noccas,C,All.hists){
  
  H <- parms$H
  firstcap <- C[H]
  pbeta <- parms$pbeta
  phibeta <- parms$phibeta
  if(DM$mod.p.h){
    zp <- parms$zp
  } else {
    zp <- rep(0,length(H))
  }
  if(DM$mod.phi.h){
    zphi <- parms$zphi
  } else {
    zphi <- rep(0,length(H))
  }
  q <- parms$q
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
  
  XBp=DM$p%*%pbeta
  XBphi=DM$phi%*%phibeta    
  
  Find <- which(firstcap<noccas)
  ind <- c(0,cumsum(seq(noccas-1,2))) 
  DMind <- rep(ind[firstcap[Find]],each=noccas-1)+unlist(lapply(firstcap[Find],function(x) c(rep(0,x-1),1:(noccas-x))))
  p <- pmin(pmax(pnorm(matrix(XBp[DMind]+rep(zp[Find],each=noccas-1),byrow=TRUE,nrow=length(Find),ncol=noccas-1)),tol),1.-tol)
  phi <- pmin(pmax(pnorm(matrix(XBphi[DMind]+rep(zphi[Find],each=noccas-1),byrow=TRUE,nrow=length(Find),ncol=noccas-1)),tol),1.-tol)
  indhists <- All.hists[H[Find],-1]
  base::sum(log( ((indhists==0) * ((1.-p) * phi *q[Find,2:noccas] + (1.-phi)*(1.-q[Find,2:noccas]))
                  + (indhists==1) * p * delta_1 * phi
                  + (indhists==2) * p * delta_2 * phi
                  + (indhists==3) * p * (1.-delta_1-delta_2) *(1.-alpha) * phi
                  + (indhists==4) * p * (1.-delta_1-delta_2) *alpha * phi)[which(q[Find,1:(noccas-1)]==1)] ))
  
}

priorsCJS<-function(parms,DM,priorparms,data_type,C,noccas){
  
  firstcap <- C[parms$H]
  
  priors <- (dmvnorm(parms$pbeta,priorparms$pbeta0,priorparms$pSigma0,log=TRUE)
             + dmvnorm(parms$phibeta,priorparms$phibeta0,priorparms$phiSigma0,log=TRUE))

  
  if(DM$mod.delta != ~NULL){
    if(DM$mod.delta==formula(~type)){
      priors <- priors + ddirichlet(c(parms$delta_1,parms$delta_2,1.-parms$delta_1-parms$delta_2),priorparms$a0delta)
    } else {
      priors <- priors + dbeta(2*parms$delta,priorparms$a0delta[1],priorparms$a0delta[2],log=TRUE)
    }
    if(data_type=="sometimes"){
      priors <- priors + dbeta(parms$alpha,priorparms$a0alpha,priorparms$b0alpha,log=TRUE)
    }
    priors <- priors + (dbeta(parms$psi,priorparms$a0psi,priorparms$b0psi,log=TRUE)
                    #+ base::sum(dbinom((firstcap<noccas),1,parms$psi,log=TRUE)))
                     + base::sum(dbinom((parms$H>1),1,parms$psi,log=TRUE)))
  }
  
  if(DM$mod.p.h){
    priors <- priors + (base::sum(dnorm(parms$zp,0.0,sqrt(parms$sigma2_zp),log=TRUE))
                        + dinvgamma(parms$sigma2_zp,shape=priorparms$l0p,scale=priorparms$d0p))
  }       
  if(DM$mod.phi.h){
    priors <- priors + (base::sum(dnorm(parms$zphi,0.0,sqrt(parms$sigma2_zphi),log=TRUE))
                        + dinvgamma(parms$sigma2_zphi,shape=priorparms$l0phi,scale=priorparms$d0phi))
  }  
  priors
}

posteriorCJS<-function(parms,DM,mms,priorparms){
  nchains<-length(parms)
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas)
  for(ichain in 1:nchains){
    temp<-parms[[ichain]]
    
    loglike <- loglikeCJS(temp,DM,noccas,mms@C,All.hists)
    
    if(!is.finite(loglike)) {
      stop(paste0("initial model likelihood is ",loglike," for chain ",ichain,". Try different initial values."))
    }
    
    posterior <- loglike + priorsCJS(temp,DM,priorparms,mms@data.type,mms@C,noccas)
    
    if(!is.finite(posterior)) {
      stop(paste("initial model posterior is",posterior,"for chain",ichain,". Try different initial values or prior parameters"))
    }
  }
}

checkCJS<-function(parms,parmlist,mms,DM,iter,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,pSigma0,phiSigma0,l0p,d0p,l0phi,d0phi,a0psi,b0psi,link){
  
  if(mms@data.type!="sometimes" & any(parms=="alpha")) stop("Parameter 'alpha' only applies to models for the 'sometimes' data type")
  
  params <- parms
  if(any(parms=="all")){
    if(mms@data.type=="sometimes"){
      params<-parmlist
    } else {
      params<-parmlist[which(parmlist!="alpha")]
    }
  } else {
    if(!all(match(params,parmlist,nomatch=0))) stop(paste0("'",params[match(params,parmlist,nomatch=0)==0],"' is not a valid parameter\n  "))
  }
  
  if(link=="logit"){
    if((bin<1 | bin>iter) & iter>0) stop("'bin' must be >0 and <",iter)
    if(taccept<=0 | taccept>1) stop ("'taccept' must be >0 and <=1")
    if(tuneadjust<=0 | tuneadjust>1) stop ("'tuneadjust' must be >0 and <=1")
  }
  
  if(thin>max(1,floor((iter-burnin+1)/2)) | thin<1) stop("'thin' must be >0 and <=",max(1,floor((iter-burnin+1)/2)))
  if(mms@ncolbasis & (maxnumbasis<1 | maxnumbasis>mms@ncolbasis)) stop("'maxnumbasis' must be between 1 and ",mms@ncolbasis)
  if(!all(c(a0delta,a0alpha,b0alpha,l0p,d0p,l0phi,d0phi,a0psi,b0psi,diag(as.matrix(pSigma0)),diag(as.matrix(phiSigma0)))>0)) stop("'a0delta', 'a0alpha', 'b0alpha', 'l0p', 'd0p', 'l0phi', 'd0phi', 'a0psi', 'b0psi', and diagonal elements of 'pSigma0' and 'phiSigma0' must be >0")
  
  
  mod.p.h<-DM$mod.p.h
  if(any(parms=="all")){
    if(!mod.p.h){
      params<-params[which(params!="zp" & params!="sigma2_zp")]
    }
  } else {
    if(!mod.p.h & (any(params=="zp") | any(params=="sigma2_zp"))) stop("Parameters 'sigma2_zp' and 'zp' only apply to individual heterogeneity models for p")
  }
  pdim<-ncol(DM$p)
  if(!pdim) stop("'mod.p' must include at least 1 parameter")
  
  mod.phi.h<-DM$mod.phi.h
  if(any(parms=="all")){
    if(!mod.phi.h){
      params<-params[which(params!="zphi" & params!="sigma2_zphi")]
    }
  } else {
    if(!mod.phi.h & (any(params=="zphi") | any(params=="sigma2_zphi"))) stop("Parameters 'sigma2_zphi' and 'zphi' only apply to individual heterogeneity models for phi")
  }
  phidim<-ncol(DM$phi)
  if(!phidim) stop("'mod.phi' must include at least 1 parameter")
  
  params
}

processCJSchains<-function(chains,params,DM,M,noccas,nchains,iter,burnin,thin){
  
  parms<-params
  if(any(parms=="phibeta")){
    parms<-c(paste0("phibeta[",colnames(DM$phi),"]"),params[which(params!="phibeta")])
  }
  if(any(parms=="pbeta")){
    parms<-c(paste0("pbeta[",colnames(DM$p),"]"),parms[which(parms!="pbeta")])
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
  if(any(parms=="zphi")){
    zphiname<-paste0("zphi[",1:M,"]")
    parms<-c(zphiname,parms[which(parms!="zphi")])
  } else {
    zphiname<-NULL
  }
  if(any(parms=="q")){
    zname<-paste0("q[",rep(1:M,each=noccas),",",1:noccas,"]")
    parms<-c(zname,parms[which(parms!="q")])
  } else {
    zname<-NULL
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
    names(initstemp) <- c(paste0("pbeta[",colnames(DM$p),"]"),paste0("phibeta[",colnames(DM$phi),"]"),"delta_1","delta_2","alpha","psi",zpname,"sigma2_zp",zphiname,"sigma2_zphi",zname,Hname,"loglike")
    if(any(params=="zp")){
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],phibeta=initstemp[paste0("phibeta[",colnames(DM$phi),"]")],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],psi=initstemp["psi"],zp=initstemp[zpname],sigma2_zp=initstemp["sigma2_zp"])
    } else {
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],phibeta=initstemp[paste0("phibeta[",colnames(DM$phi),"]")],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],psi=initstemp["psi"],zp=chains[[ichain]]$zp,sigma2_zp=initstemp["sigma2_zp"])
      names(initial.values[[ichain]]$zp) <- paste0("zp[",1:M,"]")
    }
    if(any(params=="zphi")){
      initial.values[[ichain]]$zphi <- initstemp[zphiname]
    } else {
      initial.values[[ichain]]$zphi <- chains[[ichain]]$zphi
      names(initial.values[[ichain]]$zphi) <- paste0("zphi[",1:M,"]")
    }
    initial.values[[ichain]]$sigma2_zphi <- initstemp["sigma2_zphi"]
    if(any(params=="q")){
      initial.values[[ichain]]$q <- matrix(initstemp[zname],byrow=TRUE,ncol=noccas)
    } else {
      initial.values[[ichain]]$q <- matrix(chains[[ichain]]$q,byrow=TRUE,ncol=noccas)
    }
    if(any(params=="H")){
      initial.values[[ichain]]$H <- initstemp[Hname]
    } else {
      initial.values[[ichain]]$H <- chains[[ichain]]$H
      names(initial.values[[ichain]]$H) <- paste0("H[",1:M,"]")
    }
    initial.values[[ichain]]$x <- chains[[ichain]]$x
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

#' Fit open population survival models for ``traditional'' capture-mark-recapture data consisting of a single mark type
#' 
#' This function fits Cormack-Jolly-Seber (CJS) open population models for survival probability (\eqn{\phi}) and capture probability (\eqn{p}) for ``traditional'' capture-mark-recapture data consisting of a single mark type. Using Bayesian analysis methods, Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix of observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions. With a single mark type, encounter histories consist of only non-detections (0) and type 1 encounters (1).
#' @param covs A data frame of temporal covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of sampling occasions. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mod.p Model formula for detection probability (\eqn{p}). For example, \code{mod.p=~1} specifies no effects (i.e., intercept only), \code{mod.p~time} specifies temporal effects, \code{mod.p~age} specifies age effects, \code{mod.p~h} specifies individual heterogeneity, and \code{mod.p~time+age} specifies additive temporal and age effects.
#' @param mod.phi Model formula for survival probability (\eqn{\phi}). For example, \code{mod.phi=~1} specifies no effects (i.e., intercept only), \code{mod.phi~time} specifies temporal effects, \code{mod.phi~age} specifies age effects, \code{mod.phi~h} specifies individual heterogeneity, and \code{mod.phi~time+age} specifies additive temporal and age effects.
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are probit-scale detection probability parameters ("\code{pbeta}" for \eqn{p} and "\code{phibeta}" for \eqn{\phi}), probit-scale individual heterogeneity variance terms ("\code{sigma2_zp}" for \eqn{p} and "\code{sigma2_zphi}" for \eqn{\phi}), and probit-scale individual effects ("\code{zp}" and "\code{zphi}"). Latent variable indicators for whether each individual was alive (1) or dead (0) during each sampling occasion ("\code{q}") and the log likelihood ("\code{loglike}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt Ignored; no adaptive phase is needed for "probit" link.
#' @param bin Ignored; no adaptive phase is needed for "probit" link.
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Ignored; no adaptive phase is needed for "probit" link.
#' @param tuneadjust Ignored; no adaptive phase is needed for "probit" link.
#' @param proppbeta Ignored; no adaptive phase is needed for "probit" link.
#' @param propzp Ignored; no adaptive phase is needed for "probit" link.
#' @param propsigmap Ignored; no adaptive phase is needed for "probit" link.
#' @param propphibeta Ignored; no adaptive phase is needed for "probit" link.
#' @param propzphi Ignored; no adaptive phase is needed for "probit" link.
#' @param propsigmaphi Ignored; no adaptive phase is needed for "probit" link.
#' @param pbeta0 Scaler or vector (of length k) specifying mean of pbeta ~ multivariateNormal(pbeta0, pSigma0) prior. If \code{pbeta0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{pbeta0 = 0}.  
#' @param pSigma0 Scaler or k x k matrix specifying covariance matrix of pbeta ~ multivariateNormal(pbeta0, pSigma0) prior. If \code{pSigma0} is a scaler, then this value is used for all pSigma0[j,j] for j = 1, ..., k (with pSigma[j,l] = 0 for all \eqn{j \ne l}). Default is \code{pSigma0 = 1}. 
#' @param phibeta0 Scaler or vector (of length k) specifying mean of phibeta ~ multivariateNormal(phibeta0, phiSigma0) prior. If \code{phibeta0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{phibeta0 = 0}.  
#' @param phiSigma0 Scaler or k x k matrix specifying covariance matrix of phibeta ~ multivariateNormal(phibeta0, phiSigma0) prior. If \code{phiSigma0} is a scaler, then this value is used for all phiSigma0[j,j] for j = 1, ..., k (with phiSigma[j,l] = 0 for all \eqn{j \ne l}). Default is \code{phiSigma0 = 1}. 
#' @param l0p Specifies "shape" parameter for [sigma2_zp] ~ invGamma(l0p,d0p) prior. Default is \code{l0p = 1}. 
#' @param d0p Specifies "scale" parameter for [sigma2_zp] ~ invGamma(l0p,d0p) prior. Default is \code{d0p = 0.01}. 
#' @param l0phi Specifies "shape" parameter for [sigma2_zphi] ~ invGamma(l0phi,d0phi) prior. Default is \code{l0phi = 1}. 
#' @param d0phi Specifies "scale" parameter for [sigma2_zphi] ~ invGamma(l0phi,d0phi) prior. Default is \code{d0phi = 0.01}. 
#' @param link Link function for survival and capture probabilities. Only probit link is currently implemented.
#' @param initial.values OOptional list of \code{nchain} list(s) specifying initial values for "\code{pbeta}", "\code{phibeta}", "\code{sigma2_zp}", "\code{sigma2_zphi}", "\code{zp}", "\code{zphi}", and "\code{q}". Default is \code{initial.values = NULL}, which causes initial values to be generated automatically. 
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p} and \code{mod.phi}. See \code{RMark::\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{markCJS} (or \code{\link{markClosed}}) is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{multimarkCJS}. Note that setting \code{parms="all"} is required for any \code{markCJS} model output to be used in \code{\link{multimodelCJS}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.phi}{Model formula for survival probability (as specified by \code{mod.phi} above).}
#' \item{mod.delta}{Formula always \code{NULL}; only for internal use in \code{\link{multimodelCJS}}.}
#' \item{DM}{A list of design matrices for detection and survival probability respectively generated by \code{mod.p} and \code{mod.phi}, where DM$p is the design matrix for capture probability (\eqn{p}) and DM$phi is the design matrix for survival probability (\eqn{\phi}).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{phibeta}", "\code{sigma2_zp}", "\code{sigma2_zphi}", "\code{zp}", "\code{zphi}", and "\code{q}".}
#' \item{mms}{An object of class \code{multimarksetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{processdata}}, \code{\link{multimodelCJS}}
#' @examples
#' \dontshow{
#' data<-simdataCJS(delta_1=1,delta_2=0)$Enc.Mat
#' test<-markCJS(data,iter=10,burnin=0)
#' test.age <- markCJS(data,mod.phi=~age,iter=10,burnin=0,parameters=list(Phi=list(age.bins=c(0,1,4))),right=FALSE)}
#' \donttest{
#' # These examples are excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Simulate open population data using defaults
#' data <- simdataCJS(delta_1=1,delta_2=0)$Enc.Mat
#' 
#' #Fit default open population model
#' sim.dot <- markCJS(data)
#' 
#' #Posterior summary for monitored parameters
#' summary(sim.dot$mcmc)
#' plot(sim.dot$mcmc)
#' 
#' #Fit ``age'' model with 2 age classes (e.g., juvenile and adult) for survival
#' #using 'parameters' and 'right' arguments from RMark::make.design.data
#' sim.age <- markCJS(data,mod.phi=~age,
#'            parameters=list(Phi=list(age.bins=c(0,1,4))),right=FALSE)
#' summary(getprobsCJS(sim.age))}
markCJS<-function(Enc.Mat,covs=data.frame(),mod.p=~1,mod.phi=~1,parms=c("pbeta","phibeta"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propzp=1,propsigmap=1,propphibeta=0.1,propzphi=1,propsigmaphi=1,pbeta0=0,pSigma0=1,phibeta0=0,phiSigma0=1,l0p=1,d0p=0.01,l0phi=1,d0phi=0.01,initial.values=NULL,link="probit",printlog=FALSE,...){
  if(any(Enc.Mat>1 | Enc.Mat<0)) stop("With a single mark type, encounter histories can only contain 0's (non-detections) and 1's (detections)")
  mms <- processdata(Enc.Mat,covs=covs,known=rep(1,nrow(Enc.Mat)))
  out <- multimarkCJS(mms=mms,mod.p=mod.p,mod.phi=mod.phi,mod.delta=~NULL,parms=parms,nchains=nchains,iter=iter,adapt=adapt,bin=bin,thin=thin,burnin=burnin,taccept=taccept,tuneadjust=tuneadjust,proppbeta=proppbeta,propzp=propzp,propsigmap=propsigmap,propphibeta=propphibeta,propzphi=propzphi,propsigmaphi=propsigmaphi,pbeta0=pbeta0,pSigma0=pSigma0,phibeta0=phibeta0,phiSigma0=phiSigma0,l0p=l0p,d0p=d0p,l0phi=l0phi,d0phi=d0phi,initial.values=initial.values,link=link,printlog=printlog,...)
  out$initial.values <- lapply(out$initial.values,function(x) list(pbeta=x$pbeta,phibeta=x$phibeta,zp=x$zp,sigma2_zp=x$sigma2_zp,zphi=x$zphi,sigma2_zphi=x$sigma2_zphi,q=x$q))
  return(out)
}

#' Fit open population survival models for capture-mark-recapture data consisting of multiple non-invasive marks
#' 
#' This function fits Cormack-Jolly-Seber (CJS) open population models for survival probability (\eqn{\phi}) and capture probability (\eqn{p}) from capture-mark-recapture data consisting of multiple non-invasive marks. Using Bayesian analysis methods, Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
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
#' @param mod.p Model formula for detection probability (\eqn{p}). For example, \code{mod.p=~1} specifies no effects (i.e., intercept only), \code{mod.p~time} specifies temporal effects, \code{mod.p~age} specifies age effects, \code{mod.p~h} specifies individual heterogeneity, and \code{mod.p~time+age} specifies additive temporal and age effects.
#' @param mod.phi Model formula for survival probability (\eqn{\phi}). For example, \code{mod.phi=~1} specifies no effects (i.e., intercept only), \code{mod.phi~time} specifies temporal effects, \code{mod.phi~age} specifies age effects, \code{mod.phi~h} specifies individual heterogeneity, and \code{mod.phi~time+age} specifies additive temporal and age effects.
#' @param mod.delta Model formula for conditional probabilities of type 1 (delta_1) and type 2 (delta_2) encounters, given detection. Currently only \code{mod.delta=~1} (i.e., \eqn{\delta_1 = \delta_2}) and \code{mod.delta=~type} (i.e., \eqn{\delta_1 \ne \delta_2}) are implemented.
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are probit-scale detection probability parameters ("\code{pbeta}" for \eqn{p} and "\code{phibeta}" for \eqn{\phi}), conditional probability of type 1 or type 2 encounter, given detection ("\code{delta})", probability of simultaneous type 1 and type 2 detection, given both types encountered ("\code{alpha}"), probit-scale individual heterogeneity variance terms ("\code{sigma2_zp}" for \eqn{p} and "\code{sigma2_zphi}" for \eqn{\phi}), probit-scale individual effects ("\code{zp}" and "\code{zphi}"), and the probability that a randomly selected individual from the \code{M = nrow(Enc.Mat)} observed individuals belongs to the \eqn{n} unique individuals encountered at least once ("\code{psi}"). Individual encounter history indices ("\code{H}"), latent variable indicators for whether each individual was alive (1) or dead (0) during each sampling occasion ("\code{q}"), and the log likelihood ("\code{loglike}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt Ignored; no adaptive phase is needed for "probit" link.
#' @param bin Ignored; no adaptive phase is needed for "probit" link.
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Ignored; no adaptive phase is needed for "probit" link.
#' @param tuneadjust Ignored; no adaptive phase is needed for "probit" link.
#' @param proppbeta Ignored; no adaptive phase is needed for "probit" link.
#' @param propzp Ignored; no adaptive phase is needed for "probit" link.
#' @param propsigmap Ignored; no adaptive phase is needed for "probit" link.
#' @param propphibeta Ignored; no adaptive phase is needed for "probit" link.
#' @param propzphi Ignored; no adaptive phase is needed for "probit" link.
#' @param propsigmaphi Ignored; no adaptive phase is needed for "probit" link.
#' @param maxnumbasis Maximum number of basis vectors to use when proposing latent history frequency updates. Default is \code{maxnumbasis = 1}, but higher values can potentially improve mixing.
#' @param a0delta Scaler or vector (of length d) specifying the prior for the conditional (on detection) probability of type 1 (delta_1), type 2 (delta_2), and both type 1 and type 2 encounters (1-delta_1-delta_2). If \code{a0delta} is a scaler, then this value is used for all a0delta[j] for j = 1, ..., d. For \code{mod.delta=~type}, d=3 with [delta_1, delta_2, 1-delta_1-delta_2] ~ Dirichlet(a0delta) prior. For \code{mod.delta=~1}, d=2 with [tau] ~ Beta(a0delta[1],a0delta[2]) prior, where (delta_1,delta_2,1-delta_1-delta_2) = (tau/2,tau/2,1-tau). See McClintock et al. (2013) for more details.
#' @param pbeta0 Scaler or vector (of length k) specifying mean of pbeta ~ multivariateNormal(pbeta0, pSigma0) prior. If \code{pbeta0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{pbeta0 = 0}.  
#' @param pSigma0 Scaler or k x k matrix specifying covariance matrix of pbeta ~ multivariateNormal(pbeta0, pSigma0) prior. If \code{pSigma0} is a scaler, then this value is used for all pSigma0[j,j] for j = 1, ..., k (with pSigma[j,l] = 0 for all \eqn{j \ne l}). Default is \code{pSigma0 = 1}. 
#' @param phibeta0 Scaler or vector (of length k) specifying mean of phibeta ~ multivariateNormal(phibeta0, phiSigma0) prior. If \code{phibeta0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{phibeta0 = 0}.  
#' @param phiSigma0 Scaler or k x k matrix specifying covariance matrix of phibeta ~ multivariateNormal(phibeta0, phiSigma0) prior. If \code{phiSigma0} is a scaler, then this value is used for all phiSigma0[j,j] for j = 1, ..., k (with phiSigma[j,l] = 0 for all \eqn{j \ne l}). Default is \code{phiSigma0 = 1}. 
#' @param l0p Specifies "shape" parameter for [sigma2_zp] ~ invGamma(l0p,d0p) prior. Default is \code{l0p = 1}. 
#' @param d0p Specifies "scale" parameter for [sigma2_zp] ~ invGamma(l0p,d0p) prior. Default is \code{d0p = 0.01}. 
#' @param l0phi Specifies "shape" parameter for [sigma2_zphi] ~ invGamma(l0phi,d0phi) prior. Default is \code{l0phi = 1}. 
#' @param d0phi Specifies "scale" parameter for [sigma2_zphi] ~ invGamma(l0phi,d0phi) prior. Default is \code{d0phi = 0.01}. 
#' @param a0alpha Specifies "shape1" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{a0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param b0alpha Specifies "shape2" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{b0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param a0psi Specifies "shape1" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{a0psi = 1}.
#' @param b0psi Specifies "shape2" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{b0psi = 1}.
#' @param link Link function for survival and capture probabilities. Only probit link is currently implemented.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for parameters and latent variables. Default is \code{initial.values = NULL}, which causes initial values to be generated automatically. In addition to the parameters ("\code{pbeta}", "\code{phibeta}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_zp}", "\code{sigma2_zphi}", "\code{zp}", "\code{zphi}", and "\code{psi}"), initial values can be specified for the initial latent history frequencies ("\code{x}"), initial individual encounter history indices ("\code{H}"), and initial latent variable indicators for whether each individual was alive (1) or dead (0) during each sampling occasion ("\code{q}").
#' @param known Optional integer vector indicating whether the encounter history of an individual is known with certainty (i.e., the observed encounter history is the true encounter history). Encounter histories with at least one type 4 encounter are automatically assumed to be known, and \code{known} does not need to be specified unless there exist encounter histories that do not contain a type 4 encounter that happen to be known with certainty (e.g., from independent telemetry studies). If specified, \code{known = c(v_1,v_2,...,v_M)} must be a vector of length \code{M = nrow(Enc.Mat)} where \code{v_i = 1} if the encounter history for individual \code{i} is known (\code{v_i = 0} otherwise). Note that known all-zero encounter histories (e.g., `000') are ignored.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p} and \code{mod.phi}. See \code{RMark::\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{multimarkCJS} (or \code{\link{multimarkClosed}}) is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{multimarkCJS}. Note that setting \code{parms="all"} is required for any \code{multimarkCJS} model output to be used in \code{\link{multimodelCJS}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.phi}{Model formula for survival probability (as specified by \code{mod.phi} above).}
#' \item{mod.delta}{Formula always \code{NULL}; only for internal use in \code{\link{multimodelCJS}}.}
#' \item{DM}{A list of design matrices for detection and survival probability respectively generated by \code{mod.p} and \code{mod.phi}, where DM$p is the design matrix for capture probability (\eqn{p}) and DM$phi is the design matrix for survival probability (\eqn{\phi}).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{phibeta}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_zp}" "\code{sigma2_zphi}", "\code{zp}", "\code{zphi}", "\code{psi}", "\code{x}", "\code{H}", and "\code{q}".}
#' \item{mms}{An object of class \code{multimarksetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{processdata}}, \code{\link{multimodelCJS}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#'
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#'
#' McClintock, B. T., Bailey, L. L., Dreher, B. P., and Link, W. A. 2014. Probit models for capture-recapture data subject to imperfect detection, individual heterogeneity and misidentification. \emph{The Annals of Applied Statistics} 8: 2461-2484.
#' @examples
#' \dontshow{
#' test<-multimarkCJS(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0)
#' test.age <- multimarkCJS(Enc.Mat=bobcat,mod.phi=~age,iter=10,burnin=0,parameters=list(Phi=list(age.bins=c(0,1,7))),right=FALSE)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Simulate open population data using defaults
#' data <- simdataCJS()
#' 
#' #Fit default open population model
#' sim.dot <- multimarkCJS(data$Enc.Mat)
#' 
#' #Posterior summary for monitored parameters
#' summary(sim.dot$mcmc)
#' plot(sim.dot$mcmc)
#' 
#' #' #Fit ``age'' model with 2 age classes (e.g., juvenile and adult) for survival
#' #using 'parameters' and 'right' arguments from RMark::make.design.data
#' sim.age <- multimarkCJS(data$Enc.Mat,mod.phi=~age,
#'            parameters=list(Phi=list(age.bins=c(0,1,4))),right=FALSE)
#' summary(getprobsCJS(sim.age))}
multimarkCJS<-function(Enc.Mat,data.type="never",covs=data.frame(),mms=NULL,mod.p=~1,mod.phi=~1,mod.delta=~type,parms=c("pbeta","phibeta","delta"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propzp=1,propsigmap=1,propphibeta=0.1,propzphi=1,propsigmaphi=1,maxnumbasis=1,pbeta0=0,pSigma0=1,phibeta0=0,phiSigma0=1,l0p=1,d0p=0.01,l0phi=1,d0phi=0.01,a0delta=1,a0alpha=1,b0alpha=1,a0psi=1,b0psi=1,initial.values=NULL,known=integer(),link="probit",printlog=FALSE,...){
  
  if(is.null(mms)) mms <- processdata(Enc.Mat,data.type,covs,known)
  if(class(mms)!="multimarksetup") stop("'mms' must be an object of class 'multimarksetup'")
  validObject(mms)
  
  if(class(mod.p)!="formula") stop("'mod.p' must be an object of class 'formula'")
  if(class(mod.phi)!="formula") stop("'mod.phi' must be an object of class 'formula'")
  if(class(mod.delta)!="formula") stop("'mod.delta' must be an object of class 'formula'")
  DM<-get_DMCJS(mod.p,mod.phi,mod.delta,mms@Enc.Mat,covs=mms@covs,...)
  
  if(iter>0){
    if(iter<=burnin) stop(paste("'burnin' must be less than ",iter))
  } else {
    burnin<-0
  }
  
  if(mod.delta != ~NULL) {
    parmlist<-c("pbeta","phibeta","delta","sigma2_zp","zp","sigma2_zphi","zphi","alpha","psi","q","H","loglike")
  } else {
    parmlist<-c("pbeta","phibeta","sigma2_zp","zp","sigma2_zphi","zphi","q","loglike")    
  }
  params <- checkCJS(parms,parmlist,mms,DM,iter,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,pSigma0,phiSigma0,l0p,d0p,l0phi,d0phi,a0psi,b0psi,link)
  
  data.type<-mms@data.type
  Enc.Mat<-mms@Enc.Mat
  M<-nrow(Enc.Mat)
  noccas<-ncol(Enc.Mat)
  covs<-mms@covs
  pdim<-ncol(DM$p)  
  phidim<-ncol(DM$phi)
  
  pbeta0 <-checkvecs(pbeta0,pdim,"pbeta0")
  pSigma0 <- checkmats(pSigma0,pdim,"pSigma0")  
  pprec0 <- solve(pSigma0)
  phibeta0 <- checkvecs(phibeta0,phidim,"phibeta0")
  phiSigma0 <- checkmats(phiSigma0,phidim,"phiSigma0")  
  phiprec0 <- solve(phiSigma0)
  a0delta <- checkvecs(a0delta,ifelse(mod.delta==formula(~type),3,2),"a0delta")
  
  inits<-get_inits(mms,nchains,initial.values,M,data.type,a0alpha,b0alpha,a0delta,a0psi,b0psi,DM)
  
  priorparms <-list(a0delta=a0delta,a0alpha=a0alpha,b0alpha=b0alpha,pbeta0=pbeta0,pSigma0=pSigma0,phibeta0=phibeta0,phiSigma0=phiSigma0,l0p=l0p,d0p=d0p,l0phi=l0phi,d0phi=d0phi,a0psi=a0psi,b0psi=b0psi)
  
  message("\nFitting open population model with ",link," link\n")
  if(mod.delta != ~NULL) message("data type = \"",data.type,"\"\n")
  message("p model = ",as.character(mod.p))
  message("phi model = ",as.character(mod.phi))
  if(mod.delta != ~NULL) message("delta model = ",as.character(mod.delta))
  message("\nInitializing model \n")
  posteriorCJS(inits,DM,mms,priorparms)
  
  propzp <- checkvecs(propzp,M,"propzp")
  proppbeta <- checkvecs(proppbeta,pdim,"proppbeta")
  if(length(propsigmap)!=1) stop("'propsigmap' must be a scaler")
  propzphi <- checkvecs(propzphi,M,"propzphi")
  propphibeta <- checkvecs(propphibeta,phidim,"proppbeta")
  if(length(propsigmaphi)!=1) stop("'propsigmaphi' must be a scaler")
  
  Prop.sdp <- c(propzp,proppbeta,propsigmap)
  Prop.sdphi <- c(propzphi,propphibeta,propsigmaphi)
  
  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    modlog <- ifelse(mod.delta != ~NULL,"multimarkCJS","markCJS")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0(modlog,"_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("mcmcCJS"),envir=environment())  
    chains <- parLapply(cl,1:nchains, function(ichain) mcmcCJS(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sdp,Prop.sdphi,maxnumbasis,pbeta0,pprec0,phibeta0,phiprec0,l0p,d0p,l0phi,d0phi,a0delta,a0alpha,b0alpha,a0psi,b0psi,link,printlog))
    stopCluster(cl)
    gc()
  } else {
    chains <- vector('list',nchains)
    chains[[nchains]] <- mcmcCJS(nchains,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sdp,Prop.sdphi,maxnumbasis,pbeta0,pprec0,phibeta0,phiprec0,l0p,d0p,l0phi,d0phi,a0delta,a0alpha,b0alpha,a0psi,b0psi,link,printlog)
    gc()
  }
  
  chains <- processCJSchains(chains,params,DM,M,noccas,nchains,iter,burnin,thin)
  return(list(mcmc=chains$chains,mod.p=mod.p,mod.phi=mod.phi,mod.delta=mod.delta,DM=list(p=DM$p,phi=DM$phi),initial.values=chains$initial.values,priorparms=priorparms,mms=mms))
}

#' Calculate posterior capture and survival probabilities
#'
#' This function calculates posterior capture (\eqn{p}) and survival (\eqn{\phi}) probabilities for each sampling occasion from \code{\link{multimarkCJS}} output. 
#'
#'
#' @param out List of output returned by \code{\link{multimarkCJS}}
#' @param link Link function for \eqn{p} and \eqn{\phi}. Must be "\code{probit}" or "\code{logit}". Note that \code{\link{multimarkCJS}} is currently implemented for the probit link only.
#' @return An object of class \code{\link[coda]{mcmc.list}} containing the following:
#' \item{p}{Posterior samples for capture probability (\eqn{p[c,t]}) for each release cohort (\eqn{c=1,\ldots,T-1}) and sampling occasion (\eqn{t=2,\ldots,T}).}
#' \item{phi}{Posterior samples for survival probability (\eqn{\phi[c,k]}) for each release cohort (\eqn{c=1,\ldots,T-1}) and interval (\eqn{k=1,\ldots,T-1}).}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkCJS}}
#' @examples
#' \dontshow{
#' test<-getprobsCJS(multimarkCJS(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0))}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Simulate open population data with temporal variation in survival
#' noccas <- 5
#' data <- simdataCJS(noccas=noccas, phibeta=rnorm(noccas-1,1.6,0.1))
#'  
#' #Fit open population model with temporal variation in survival
#' sim.time <- multimarkCJS(data$Enc.Mat,mod.phi=~time)
#'     
#' #Calculate capture and survival probabilities for each cohort and time
#' pphi <- getprobsCJS(sim.time)
#' summary(pphi)}
getprobsCJS<-function(out,link="probit"){
  
  DMp<-out$DM$p
  DMphi<-out$DM$phi
  mod.p.h<-any("h"==attributes(terms(out$mod.p))$term.labels)
  mod.phi.h<-any("h"==attributes(terms(out$mod.phi))$term.labels)
  
  noccas<-ncol(out$initial.values[[1]]$q)
  if(noccas<2) stop("must have >1 sampling occasion")
  
  pbetanames<-paste0("pbeta[",colnames(DMp),"]")
  phibetanames<-paste0("phibeta[",colnames(DMphi),"]")
  nchains<-length(out$mcmc)
  
  pphi<-vector("list",nchains)
  
  varind <- is.null(varnames(out$mcmc))
  if(!varind){
    vars <- varnames(out$mcmc)
  } else {
    vars <- names(out$mcmc[[1]])    
  }
  if(!any(match(pbetanames,vars,nomatch=0))) stop("'pbeta' parameters not found")
  if(!any(match(phibetanames,vars,nomatch=0))) stop("'phibeta' parameters not found")
  
  for(ichain in 1:nchains){
    
    p <- inverseXB(ichain,out,pbetanames,mod.p.h,DMp,nrow(DMp),varind,vars,"p","sigma2_zp",link)
    phi <- inverseXB(ichain,out,phibetanames,mod.phi.h,DMphi,nrow(DMphi),varind,vars,"phi","sigma2_zphi",link)
    
    pphiseq<-matrix(0,(noccas-1),(noccas-1))
    for(t in 1:(noccas-1)){
      pphiseq[t:(noccas-1),t]<-t:(noccas-1)
    }
    colnames(p) <- paste0("p[",rep(1:(noccas-1),times=(noccas-1):1),",",pphiseq[which(pphiseq>0)]+1,"]")
    colnames(phi) <- paste0("phi[",rep(1:(noccas-1),times=(noccas-1):1),",",pphiseq[which(pphiseq>0)],"]")
    
    pphi[[ichain]]<- mcmc(cbind(p,phi),start=start(out$mcmc),end=end(out$mcmc),thin=attributes(out$mcmc[[ichain]])$mcpar[3])
  }
  return(as.mcmc.list(pphi))
}

checkparmsCJS <- function(mms,modlist,params,parmlist,M){    
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
    if((length(deltatypeind)+length(delta1ind))!=length(modlist)) stop("Cannot perform multimodel inference using both 'multimarkCJS()' and 'markCJS()' models")
  }
  hpind <- which(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.p))$term.labels))==TRUE)  
  if(!length(hpind)){
    if(!all(lapply(params,function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
  } else {
    if(!all(lapply(params[-hpind],function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
    tmpparmlist<-c(parmlist,"sigma2_zp",paste0("zp[",1:M,"]"))
    if(!all(lapply(params[hpind],function(x) base::sum(match(x,tmpparmlist,nomatch=0)))==base::sum(1:length(tmpparmlist))))  stop("required parameters not found for all models")
  }
  hphiind <- which(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.phi))$term.labels))==TRUE)  
  if(!length(hphiind)){
    if(!all(lapply(params,function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
  } else {
    if(!all(lapply(params[-hphiind],function(x) base::sum(match(x,parmlist,nomatch=0)))==base::sum(1:length(parmlist)))) stop("required parameters not found for all models")
    tmpparmlist<-c(parmlist,"sigma2_zphi",paste0("zphi[",1:M,"]"))
    if(!all(lapply(params[hphiind],function(x) base::sum(match(x,tmpparmlist,nomatch=0)))==base::sum(1:length(tmpparmlist))))  stop("required parameters not found for all models")
  }
}

checkmmCJSinput<-function(mmslist,modlist,nmod,nchains,iter,miter,mburnin,mthin,modprior,M1){
  if(!all(match(unlist(unique(lapply(modlist,names))),c("mcmc","mod.p","mod.phi","mod.delta","DM","initial.values","priorparms","mms"),nomatch=0))) stop("each object in 'modlist' must be a list returned by multimarkCJS()")
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

drawmissingCJS<-function(M.cur,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale){
  missingpbeta <- rnorm(length(missing$missingpbetaparms[[M.cur]]),sd=pbetapropsd)
  names(missingpbeta) <- missing$missingpbetaparms[[M.cur]]
  missingphibeta <- rnorm(length(missing$missingphibetaparms[[M.cur]]),sd=phibetapropsd)
  names(missingphibeta) <- missing$missingphibetaparms[[M.cur]]
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
  missingsigphi  <- rinvgamma(length(missing$missingsigphiparms[[M.cur]]),shape=sigphipropshape,scale=sigphipropscale)
  names(missingsigphi) <- missing$missingsigphiparms[[M.cur]]
  missingzp <- rnorm(length(missing$missingzpparms[[M.cur]]),sd=missing$zppropsd+sqrt(missingsigp)*missing$usesigp)
  names(missingzp) <- missing$missingzpparms[[M.cur]]
  missingzphi <- rnorm(length(missing$missingzphiparms[[M.cur]]),sd=missing$zphipropsd+sqrt(missingsigphi)*missing$usesigphi)
  names(missingzphi) <- missing$missingzphiparms[[M.cur]]
  missing <- c(missingpbeta,missingphibeta,missingdelta,missingzp,missingzphi,missingsigp,missingsigphi)
  missing
}

getbrobprobCJS<-function(imod,modprior,posterior,cur.parms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale){
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
       + base::sum(dnorm(cur.parms[missing$missingphibetaparms[[imod]]],sd=phibetapropsd,log=TRUE))
       + base::sum(dnorm(cur.parms[missing$missingzpparms[[imod]]],sd=missing$zppropsd+sqrt(cur.parms[missing$missingsigpparms[[imod]]])*missing$usesigp,log=TRUE))
       + base::sum(dnorm(cur.parms[missing$missingzphiparms[[imod]]],sd=missing$zphipropsd+sqrt(cur.parms[missing$missingsigphiparms[[imod]]])*missing$usesigphi,log=TRUE))
       + base::sum(dinvgamma(cur.parms[missing$missingsigpparms[[imod]]],shape=sigppropshape,scale=sigppropscale))
       + base::sum(dinvgamma(cur.parms[missing$missingsigphiparms[[imod]]],shape=sigphipropshape,scale=sigphipropscale))
       + deltadens)
}

getcurCJSparmslist<-function(cur.parms,DM,M,noccas,data_type,alpha){
  
  parmslist=vector('list',1)
  parmslist[[1]]$H<-cur.parms[paste0("H[",1:M,"]")]
  parmslist[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
  parmslist[[1]]$phibeta <- cur.parms[paste0("phibeta[",colnames(DM$phi),"]")]
  parmslist[[1]]$zp <- cur.parms[paste0("zp[",1:M,"]")]
  parmslist[[1]]$sigma2_zp <- cur.parms["sigma2_zp"]
  parmslist[[1]]$zphi <- cur.parms[paste0("zphi[",1:M,"]")]
  parmslist[[1]]$sigma2_zphi <- cur.parms["sigma2_zphi"]
  
  parmslist[[1]]$psi <- cur.parms["psi"]
  parmslist[[1]]$delta_1 <- cur.parms["delta_1"]
  parmslist[[1]]$delta_2 <- cur.parms["delta_2"]
  parmslist[[1]]$delta <- cur.parms["delta"]
  
  if(data_type=="sometimes"){
    parmslist[[1]]$alpha <- cur.parms["alpha"]
  } else {
    parmslist[[1]]$alpha <- alpha   
  }
  parmslist[[1]]$q <- matrix(cur.parms[paste0("q[",rep(1:M,each=noccas),",",1:noccas,"]")],byrow=T,ncol=noccas) 
  parmslist
}

missingparmnamesCJS<-function(params,M,noccas,zppropsd,zphipropsd){
  
  multiparms <- unique(unlist(params))
  
  commonparms <- Reduce(intersect, params)
  commonparms <- commonparms[-match(c(paste0("q[",rep(1:M,each=noccas),",",1:noccas,"]"),"loglike"),commonparms)]
  
  missingparms <- lapply(params,get_missingparms,multiparms=multiparms)
  
  missingpbetaparms <- extractmissingparms(missingparms,"pbeta")
  
  missingsigpparms <- lapply(missingparms,function(x) unlist(x,use.names=FALSE)[which(x=="sigma2_zp")])
  
  missingzpparms <- extractmissingparms(missingparms,"zp[")
  
  missingphibetaparms <- extractmissingparms(missingparms,"phibeta")
  
  missingsigphiparms <- lapply(missingparms,function(x) unlist(x,use.names=FALSE)[which(x=="sigma2_zphi")])
  
  missingzphiparms <- extractmissingparms(missingparms,"zphi[")
  
  missingdeltaparms <- extractmissingparms(missingparms,"delta")
  
  if(is.null(zppropsd)){
    zppropsd <- 0
    usesigp <- 1
  } else {
    usesigp <-0
  }
  
  if(is.null(zphipropsd)){
    zphipropsd <- 0
    usesigphi <- 1
  } else {
    usesigphi <-0
  }
  list(commonparms=commonparms,missingparms=missingparms,missingpbetaparms=missingpbetaparms,missingsigpparms=missingsigpparms,missingzpparms=missingzpparms,missingphibetaparms=missingphibetaparms,missingsigphiparms=missingsigphiparms,missingzphiparms=missingzphiparms,missingdeltaparms=missingdeltaparms,zppropsd=zppropsd,usesigp=usesigp,zphipropsd=zphipropsd,usesigphi=usesigphi) 
}

monitorparmsCJS <- function(parms,parmlist,noccas){ 
  
  if(!all(match(parms,parmlist,nomatch=0))) stop(paste0("monitored parameters ('monparms') can only include: ",paste(parmlist[-length(parmlist)],collapse=", "),", or ",parmlist[length(parmlist)]))
  
  commonparms <- parms
  
  getprobitp <- derivedprobitfun(parms,"p")
  getprobitphi <- derivedprobitfun(parms,"phi")   
  
  pphiseq<-matrix(0,(noccas-1),(noccas-1))
  for(t in 1:(noccas-1)){
    pphiseq[t:(noccas-1),t]<-t:(noccas-1)
  }
  
  if(any(parms=="p")){
    namesp <- paste0("p[",rep(1:(noccas-1),times=(noccas-1):1),",",pphiseq[which(pphiseq>0)]+1,"]")
    commonparms <- commonparms[-which(parms=="p")]
    parms <- parms[-which(parms=="p")]
    parms <- c(parms,namesp)
  } else {
    namesp <- NULL
  }
  if(any(parms=="phi")){
    namesphi <- paste0("phi[",rep(1:(noccas-1),times=(noccas-1):1),",",pphiseq[which(pphiseq>0)],"]")
    commonparms <- commonparms[-which(parms=="phi")]
    parms <- parms[-which(parms=="phi")]
    parms <- c(parms,namesphi)
  } else {
    namesphi <- NULL
  }
  list(commonparms=commonparms,parms=parms,namesp=namesp,namesphi=namesphi,getprobitp=getprobitp,getprobitphi=getprobitphi)
}

rjmcmcCJS <- function(ichain,mms,M,noccas,data_type,alpha,C,All.hists,modlist,DMlist,deltalist,priorlist,mod.p.h,mod.phi.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale,pmodnames,phimodnames,deltamodnames,printlog){
  
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
  
  modmissingparms <- drawmissingCJS(M.cur,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale)
  cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
  
  DM <- DMlist[[M.cur]]
  DM$mod.delta <- deltalist[[M.cur]]
  DM$mod.p.h <- mod.p.h[[M.cur]]
  DM$mod.phi.h <- mod.phi.h[[M.cur]]
  
  cur.parms.list <- getcurCJSparmslist(cur.parms,DM,M,noccas,data_type,alpha)  
  
  for(iiter in 1:miter){
    
    posterior <- cur.parms["loglike"] + priorsCJS(cur.parms.list[[1]],DM,priorlist[[M.cur]],data_type,C,noccas)
    
    mod.prob.brob[M.cur] <- getbrobprobCJS(M.cur,modprior,posterior,cur.parms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale)
    
    for(imod in (1:nmod)[-M.cur]){ 
      
      DM <- DMlist[[imod]]
      DM$mod.delta <- deltalist[[imod]]
      DM$mod.p.h <- mod.p.h[imod]
      DM$mod.phi.h <- mod.phi.h[imod]
      
      cur.parms.list[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
      cur.parms.list[[1]]$phibeta <- cur.parms[paste0("phibeta[",colnames(DM$phi),"]")]
      
      loglike <- loglikeCJS(cur.parms.list[[1]],DM,noccas,C,All.hists)
      
      posterior <- loglike + priorsCJS(cur.parms.list[[1]],DM,priorlist[[imod]],data_type,C,noccas)
      
      mod.prob.brob[imod] <- getbrobprobCJS(imod,modprior,posterior,cur.parms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale)
    }
    
    if(any(is.na(as.numeric(mod.prob.brob)))){
      warning(paste0("'NA' posterior for model '","p(",pmodnames[is.na(as.numeric(mod.prob.brob))],")phi(",phimodnames[is.na(as.numeric(mod.prob.brob))],")delta(",deltamodnames[is.na(as.numeric(mod.prob.brob))],")' at iteration ",iiter,"; model move rejected."))
      flush.console()
    } else {       
      mod.prob <- as.numeric(mod.prob.brob/Brobdingnag::sum(mod.prob.brob))
      M.cur <- (1:nmod)[rmultinom(1, 1, mod.prob)==1]
    }
    
    modmissingparms <- drawmissingCJS(M.cur,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale)
    cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
    
    DM <- DMlist[[M.cur]]
    DM$mod.delta <- deltalist[[M.cur]]
    DM$mod.p.h <- mod.p.h[[M.cur]]
    DM$mod.phi.h <- mod.phi.h[[M.cur]]
    
    cur.parms.list <- getcurCJSparmslist(cur.parms,DM,M,noccas,data_type,alpha)  
    
    if(iiter>mburnin & !iiter%%mthin){
      multimodel[iiter/mthin-floor(mburnin/mthin),"M"] <- M.cur
      multimodel[iiter/mthin-floor(mburnin/mthin),commonparms] <- cur.parms[commonparms]
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesp] <- monitorparms$getprobitp(DM$mod.p.h,DM$p,cur.parms.list[[1]]$pbeta,cur.parms.list[[1]]$sigma2_zp)
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesphi] <- monitorparms$getprobitphi(DM$mod.phi.h,DM$phi,cur.parms.list[[1]]$phibeta,cur.parms.list[[1]]$sigma2_zphi)
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

#' Multimodel inference for 'multimark' open population survival models
#' 
#' This function performs Bayesian multimodel inference for a set of 'multimark' open population survival (i.e., Cormack-Jolly-Seber) models using the reversible jump Markov chain Monte Carlo (RJMCMC) algorithm proposed by Barker & Link (2013).
#'
#'
#' @param modlist A list of individual model output lists returned by \code{\link{multimarkCJS}}. The models must have the same number of chains and MCMC iterations.
#' @param modprior Vector of length \code{length(modlist)} containing prior model probabilities. Default is \code{modprior = rep(1/length(modlist), length(modlist))}.
#' @param monparms Parameters to monitor. Only parameters common to all models can be monitored (e.g., "\code{pbeta[(Intercept)]}", "\code{phibeta[(Intercept)]}", "\code{psi}"), but derived survival ("\code{phi}") and capture ("\code{p}") probabilities can also be monitored. Default is \code{monparms = "phi"}.
#' @param miter The number of RJMCMC iterations per chain. If \code{NULL}, then the number of MCMC iterations for each individual model chain is used.
#' @param mburnin Number of burn-in iterations (\code{0 <= mburnin < miter}).
#' @param mthin Thinning interval for monitored parameters.
#' @param M1 Integer vector indicating the initial model for each chain, where \code{M1_j=i} initializes the RJMCMC algorithm for chain j in the model corresponding to \code{modlist[[i]]} for i=1,...,  \code{length(modlist)}. If \code{NULL}, the algorithm for all chains is initialized in the most general model. Default is \code{M1=NULL}.
#' @param pbetapropsd Scaler specifying the standard deviation of the Normal(0, pbetapropsd) proposal distribution for "\code{pbeta}"  parameters. Default is \code{pbetapropsd=1}. See Barker & Link (2013) for more details.
#' @param zppropsd Scaler specifying the standard deviation of the Normal(0, zppropsd) proposal distribution for "\code{zp}"  parameters. Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. If \code{NULL}, zppropsd = sqrt(sigma2_zp) is used. Default is \code{zppropsd=NULL}. See Barker & Link (2013) for more details.  
#' @param phibetapropsd Scaler specifying the standard deviation of the Normal(0, phibetapropsd) proposal distribution for "\code{phibeta}"  parameters. Default is \code{phibetapropsd=1}. See Barker & Link (2013) for more details.
#' @param zphipropsd Scaler specifying the standard deviation of the Normal(0, zphipropsd) proposal distribution for "\code{zphi}"  parameters. Only applies if at least one (but not all) model(s) include individual hetergeneity in survival probability. If \code{NULL}, zphipropsd = sqrt(sigma2_zphi) is used. Default is \code{zphipropsd=NULL}. See Barker & Link (2013) for more details.  
#' @param sigppropshape Scaler specifying the shape parameter of the invGamma(shape = sigppropshape, scale = sigppropscale) proposal distribution for "\code{sigma2_zp}". Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. Default is \code{sigppropshape=1}. See Barker & Link (2013) for more details.
#' @param sigppropscale Scaler specifying the scale parameter of the invGamma(shape = sigppropshape, scale = sigppropscale) proposal distribution for "\code{sigma2_zp}". Only applies if at least one (but not all) model(s) include individual hetergeneity in detection probability. Default is \code{sigppropscale=0.01}. See Barker & Link (2013) for more details.
#' @param sigphipropshape Scaler specifying the shape parameter of the invGamma(shape = sigphipropshape, scale = sigphipropscale) proposal distribution for "\code{sigma2_zphi}". Only applies if at least one (but not all) model(s) include individual hetergeneity in survival probability. Default is \code{sigphipropshape=1}. See Barker & Link (2013) for more details.
#' @param sigphipropscale Scaler specifying the scale parameter of the invGamma(shape = sigphipropshape, scale = sigphipropscale) proposal distribution for "\code{sigma_zphi}". Only applies if at least one (but not all) model(s) include individual hetergeneity in survival probability. Default is \code{sigphipropscale=0.01}. See Barker & Link (2013) for more details.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @details Note that setting \code{parms="all"} is required when fitting individual \code{\link{multimarkCJS}} models to be included in \code{modlist}.
#' @return A list containing the following:
#' \item{rjmcmc}{Reversible jump Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}. Includes RJMCMC output for monitored parameters and the current model at each iteration ("\code{M}").}
#' \item{pos.prob}{A list of calculated posterior model probabilities for each chain, including the overall posterior model probabilities across all chains.}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkCJS}}, \code{\link{processdata}}
#' @references
#' Barker, R. J. and Link. W. A. 2013. Bayesian multimodel inference by RJMCMC: a Gibbs sampling approach. The American Statistician 67: 150-156.
#' @examples
#' \dontshow{
#' setup<-processdata(bobcat)
#' test.dot<-multimarkCJS(mms=setup,parms="all",iter=10,burnin=0)
#' test<-multimodelCJS(modlist=list(mod1=test.dot,mod2=test.dot))
#' set.seed(10)
#' }
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Generate object of class "multimarksetup" from simulated data
#' data_type = "always"
#' noccas <- 5
#' phibetaTime <- seq(2,0,length=noccas-1) # declining trend in survival
#' data <- simdataCJS(noccas=5,phibeta=phibetaTime,data.type=data_type)
#' setup <- processdata(data$Enc.Mat,data.type=data_type)
#' 
#' #Run single chain using the default model. Note parms="all".
#' sim.pdot.phidot <- multimarkCJS(mms=setup,parms="all")
#' 
#' #Run single chain with temporal trend for phi. Note parms="all".
#' sim.pdot.phiTime <- multimarkCJS(mms=setup,mod.phi=~Time,parms="all")
#' 
#' #Perform RJMCMC using defaults
#' modlist <- list(mod1=sim.pdot.phidot,mod2=sim.pdot.phiTime)
#' sim.M <- multimodelCJS(modlist=modlist)
#' 
#' #Posterior model probabilities
#' sim.M$pos.prob
#' 
#' #multimodel posterior summary for survival (display first cohort only)
#' summary(sim.M$rjmcmc[,paste0("phi[1,",1:(noccas-1),"]")])}
multimodelCJS<-function(modlist,modprior=rep(1/length(modlist),length(modlist)),monparms="phi",miter=NULL,mburnin=0,mthin=1,M1=NULL,pbetapropsd=1,zppropsd=NULL,phibetapropsd=1,zphipropsd=NULL,sigppropshape=1,sigppropscale=0.01,sigphipropshape=1,sigphipropscale=0.01,printlog=FALSE){
  
  nmod <- length(modlist)
  iter <- unlist(unique(lapply(modlist,function(x) unique(lapply(x$mcmc,nrow)))))
  nchains <- unlist(unique(lapply(modlist,function(x) length(x$mcmc))))
  mmslist <- unlist(unique(lapply(modlist, function(x) {x$mms@covs<-data.frame();x$mms})))
  
  params <- lapply(modlist,function(x) varnames(x$mcmc))
  
  if(is.null(M1)) M1 <- rep(which.max(lapply(params,length))[1],nchains)
  
  if(is.null(miter)) miter <- iter
  
  mms<-checkmmCJSinput(mmslist,modlist,nmod,nchains,iter,miter,mburnin,mthin,modprior,M1)
  
  noccas<-ncol(mms@Enc.Mat)
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas)
  C<-mms@C
  
  checkparmsCJS(mms,modlist,params,parmlist=c("pbeta[(Intercept)]","phibeta[(Intercept)]",paste0("q[",rep(1:M,each=noccas),",",1:noccas,"]"),"loglike"),M)
  
  pmodnames <- unlist(lapply(modlist,function(x) x$mod.p)) 
  phimodnames <- unlist(lapply(modlist,function(x) x$mod.phi))
  deltamodnames <- unlist(lapply(modlist,function(x) x$mod.delta)) 
  
  message("\nPerforming open population Bayesian multimodel inference by RJMCMC \n")
  if(all(deltamodnames!=~NULL)) {
    message(paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")delta(",deltamodnames,")\n"))  
  } else if(all(deltamodnames==~NULL)){
    message(paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")\n"))
  }
  
  missing <- missingparmnamesCJS(params,M,noccas,zppropsd,zphipropsd)
  
  monitorparms <- monitorparmsCJS(monparms,c(missing$commonparms,"p","phi"),noccas)
    
  DMlist <- lapply(modlist,function(x) x$DM)
  deltalist <- lapply(modlist,function(x) x$mod.delta)
  priorlist <- lapply(modlist,function(x) x$priorparms)
  mod.p.h <- unlist(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.p))$term.labels)))
  mod.phi.h <- unlist(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.phi))$term.labels)))
  
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
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0("multimodelCJS_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("rjmcmcCJS"),envir=environment())  
    multimodel <- parLapply(cl,1:nchains, function(ichain) rjmcmcCJS(ichain,mms,M,noccas,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[ichain]]),DMlist,deltalist,priorlist,mod.p.h,mod.phi.h,iter,miter,mburnin,mthin,modprior,M1[ichain],monitorparms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale,pmodnames,phimodnames,deltamodnames,printlog))
    stopCluster(cl)
    gc()
  } else {
    multimodel <- vector('list',nchains)
    multimodel[[nchains]] <- rjmcmcCJS(nchains,mms,M,noccas,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[nchains]]),DMlist,deltalist,priorlist,mod.p.h,mod.phi.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,phibetapropsd,sigppropshape,sigppropscale,sigphipropshape,sigphipropscale,pmodnames,phimodnames,deltamodnames,printlog)
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
      names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")delta(",deltamodnames,")")
    } else {
      names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")")
    }
    multimodel[[ichain]] <- mcmc(multimodel[[ichain]])
    attributes(multimodel[[ichain]])$mcpar <- c(head(temp,n=1),tail(temp,n=1),mthin)
  }
  
  multimodel <- as.mcmc.list(multimodel)
  names(pos.prob) <- paste0("chain",1:nchains)
  pos.prob[["overall"]]<- hist(unlist(multimodel[, "M"]),plot = F, breaks = 0:nmod)$density
  if(all(deltamodnames!=~NULL)){
    names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")delta(",deltamodnames,")")
  } else {
    names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")phi(",phimodnames,")")
  }
  list(rjmcmc=multimodel,pos.prob=pos.prob) 
}