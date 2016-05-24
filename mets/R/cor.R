dep.cif<-function(cif,data,cause=NULL,model="OR",cif2=NULL,times=NULL,
                  cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
                  clusters=NULL,theta=NULL,theta.des=NULL,step=1,sym=1,weights=NULL,
		  same.cens=FALSE,censoring.weights=NULL,silent=1,entry=NULL,estimator=1,
		  trunkp=1,admin.cens=NULL,control=list(),par.func=NULL,dpar.func=NULL,dimpar=NULL,
		  score.method="nlminb",random.design=NULL,exp.link=0,...)
{ ## {{{
  ## {{{ set up data and design
  multi<-0; dscore=1; stab.cens<-FALSE; entry.call<-entry; inverse<-exp.link
  notaylor<-1; flex.func<-0; 

  ## extract design and time and cause from cif object 
  time <- cif$response[,"exit"]
  if (is.null(cause)) cause <- cif$response[,"cause"]  ##  attr(cif,"cause"); 
  cause <- as.numeric(cause)
  if (is.null(cens.code)) cens.code <- attr(cif,"cens.code")
  delta<-(cause!=cens.code)
  if (length(cause)!=length(time)) stop("cause and time not of same length\n"); 
  formula <- attr(cif,"Formula")
  ldata <- aalen.des2(formula(delete.response(terms(formula))),data=data,model="aalen")
  X<-ldata$X;  Z<-ldata$Z;  

  antpers<-nrow(X); 
  if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else {Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
  if (npar==TRUE) {Z<-matrix(0,antpers,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-ncol(Z);} 
  ng<-antpers;px<-ncol(X);  
  if (is.null(times)) times<-cif$cum[,1]; 
  est<-as.matrix(cif$cum); 
  if (semi==1) gamma<-cif$gamma  else gamma<-0; 
  ntimes<-length(times); 
  if ((cif$model!="additive") && (cif$model!="fg")) 
    stop("Marginal Cumulative incidence model, must be either additive or extended Fine-Gray model\n")
  cif.model <- switch(cif$model,additive=1,fg=2)

  if (cause1!=attr(cif,"causeS")) cat("Cause for marginal model and correlation not the same\n"); 
  if ((cause1[1]!=cause2[1])) {
    if (is.null(cif2)==TRUE) stop("Must provide marginal model for both causes"); 
    formula2<-attr(cif2,"Formula"); 
    ldata2 <- aalen.des2(formula(delete.response(terms(formula2))),data=data,model="aalen");
    X2<-ldata2$X; Z2<-ldata$Z;  
    if (is.null(Z2)==TRUE) {npar2<-TRUE; semi2<-0;}  else {Z2<-as.matrix(Z2); npar2<-FALSE; semi2<-1;}
    if (npar2==TRUE) {Z2<-matrix(0,antpers,1); pg2<-1; fixed2<-0;} else {fixed2<-1;pg2<-ncol(Z2);} 
    px2<-ncol(X2);  
    est2<-cif2$cum; 
    if (semi2==1) gamma2<-cif2$gamma  else gamma2<-0; 
    est2<-Cpred(est2,times,strict=FALSE);
  } else { 
    X2<-matrix(0,1,1); Z2<-matrix(0,1,1); pg2<-1; px2<-1;  semi2<-0; est2<-matrix(0,1,2); gamma2<-0; 
    npar2<-FALSE; 
  }


### For truncation 
  if (is.null(entry)) entry.call <- NULL else entry.call <- 0
  if (is.null(entry)) entry <- rep(0,antpers); 
  cum1<-Cpred(rbind(rep(0,px+1),cif$cum),entry)[,-1];
  if (cif.model==1) cif1entry  <-  1-exp(-apply(X*cum1,1,sum)- (Z %*% gamma )*entry)
  else if (cif.model==2) cif1entry  <-  1-exp(-apply(X*cum1,1,sum)*exp(Z %*% gamma ))
  if ((model!="COR") && (cause1[1]!=cause2[1])) {
    cum2<-Cpred(rbind(rep(0,px2+1),cif2$cum),entry)[,-1];
    if (cif.model==1) cif2entry  <-  1-exp(-apply(X2*cum2,1,sum)- (Z2 %*% gamma2 )*entry)
    else if (cif.model==2) cif2entry  <-  1-exp(-apply(X2*cum2,1,sum)*exp(Z2 %*% gamma2 ))
  } else {cum2 <- cum1; cif2entry <- cif1entry;}
  ## }}}


  ## {{{ censoring model stuff
  cens.weight <- cif$cens.weights ### censoring weights from cif function
  Gcxe <- 1; 
  ### when censoring probs given, wants them so cens.model="user.weights"
  if (!is.null(censoring.weights))  cens.model <- "user.weights"
  if (cens.model!="user.weights") {
    if (is.null(cens.weight) ) {
      if (cens.model=="KM") { ## {{{
        ud.cens<-survival::survfit(Surv(time,cause==cens.code)~+1); 
        Gfit<-cbind(ud.cens$time,ud.cens$surv)
        Gfit<-rbind(c(0,1),Gfit); 
        Gcx<-Cpred(Gfit,time)[,2];
        if (is.null(entry.call)==FALSE) Gcxe<-Cpred(Gfit,entry)[,2];
        Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
        Gctimes<- Gcx ## }}}
      } else if (cens.model=="cox") { ## {{{
        if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
        ud.cens<-cox.aalen(Surv(time,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
        Gcx<-Cpred(ud.cens$cum,time)[,2];
        RR<-exp(XZ %*% ud.cens$gamma)
        Gcx<-exp(-Gcx*RR)
        Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
        if (is.null(entry.call)==FALSE) {
          Gcxe<-Cpred(Gfit,entry)[,2];
          Gcxe<-exp(-Gcxe*RR)
        }
        Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
        Gctimes<- Gcx  ## }}}
      } else if (cens.model=="aalen") {  ## {{{
        if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
        ud.cens<-aalen(Surv(time,cause==cens.code)~XZ,n.sim=0,robust=0);
        Gcx<-Cpred(ud.cens$cum,time)[,-1];
        XZ<-cbind(1,XZ); 
        Gcx<-exp(-apply(Gcx*XZ,1,sum))
        Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
        Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
        if (is.null(entry.call)==FALSE) {
          Gcxe<-Cpred(Gfit,entry)[,-1];
          Gcxe<-exp(-apply(Gcxe*XZ,1,sum))
          Gcxe[Gcxe>1]<-1; Gcxe[Gcxe<0]<-0
        }
        Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
        Gctimes<- Gcx  ## }}}
      }
    } else { Gcx <- cens.weight; Gctimes <- cens.weight;} 
  } else {
    Gcx <- censoring.weights
    Gctimes <- cens.weight
  } 
  ntimes<-length(times); 
  ## }}}

  ## {{{ set up cluster + theta design + define iid variables 
  if (is.null(clusters)== TRUE) {
    ## take clusters from cif model 
    clusters <- attr(cif,"clusters"); 
    antclust<- length(unique(clusters)); 
    max.clust <- attr(cif,"max.clust")
    if (attr(cif,"coarse.clust")) 
      stop("Max.clust should be NULL in marginal model, or clusters should be given at call \n"); 
  } else {
    clus<-unique(clusters); antclust<-length(clus);
    clusters <- as.integer(factor(clusters, labels = 1:(antclust)))-1;
  }

  outc <- cluster.index(clusters,index.type=TRUE); 
  clustsize <- outc$cluster.size
  maxclust <- outc$maxclust
  clusterindex <- outc$idclust
  if (maxclust==1) stop("No clusters given \n"); 

 if (model!="ARANCIF") { 
  if (is.null(theta.des)==TRUE) { ptheta<-1; theta.des<-matrix(1,antpers,ptheta);
  colnames(theta.des) <- "intercept"; } else theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des);
  if (is.null(dimpar) & is.null(par.func)) dimpar<-ptheta; 
 }

 if (model=="ARANCIF") { ### different parameters for Additive random effects 
     if (is.null(random.design)) random.design <- matrix(1,antpers,1); 
     dim.rv <- ncol(random.design); 
     if (is.null(theta.des)==TRUE) theta.des<-diag(dim.rv);
     dimpar <- ncol(theta.des); 
 
     if (nrow(theta.des)!=ncol(random.design)) stop("nrow(theta.des)!= ncol(random.design)"); 
###     score.method <- "nlminb"; ### force nlminb because derivatives are not working
 } else random.design <- matrix(0,1,1); 

  if ( (!is.null(par.func)) && is.null(dimpar) ) 
    stop("Must specify dimension of parameters when specifying R-functions\n")

  if ( (is.null(theta)==TRUE) && (model=="ARANCIF")) theta<-rep(0.5,dimpar);
  if (is.null(theta)==TRUE) theta<-rep(0.1,dimpar);
  if (length(theta)!=dimpar) theta<-rep(theta[1],dimpar);

  Biid<-c(); gamma.iid <- 0; B2iid<-c(); gamma2.iid <- 0; 
  if (notaylor==0) {
    for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
    if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
    if ((model!="COR") && (cause1!=cause2))  {
      B2iid<-c()
      for (i in 1:antclust) B2iid<-cbind(B2iid,cif2$B.iid[[i]]); 
      if (npar2==TRUE) gamma2.iid<-0 else  gamma2.iid<-cif2$gamma.iid;
    } else { B2iid<-Biid; gamma2.iid<-gamma.iid; }
  }
  var.theta<-hess<-matrix(0,dimpar,dimpar); score<-rep(0,dimpar); 
  time.pow<-attr(cif,"time.pow"); 
  #if (sum(time.pow)==0) time.pow<-rep(1,pg); 
  theta.iid<-matrix(0,antclust,dimpar); 

  if ((nrow(theta.des)!=nrow(X)) && (model!="ARANCIF")) 
	  stop("Dependence design not consistent with data\n"); 
  if (length(trunkp)!=antpers) trunkp <- rep(1,antpers)
  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
  ## }}}

###  ## {{{ function and derivative
  if (is.null(par.func)==FALSE)
  {
     flex.func<-1; # use flexible design
  if (is.null(dpar.func)) {
	  cat("Must provide derivative that is used for iid decomposition for SE's\n"); 
	  score.method <- "nlminb"
	  dpar.func <- par.func
  }
###  if (score.method=="fisher.scoring") 
###  cat("Score.method set to nlminb for flexible modelling \n"); 
  } ## }}}

  Zgamma <-  c(Z %*% gamma); 
  Z2gamma2 <- c(Z2 %*% gamma2); 
  dep.model <- 0;  
  dep.model <- switch(model,COR=1,RR=2,OR=3,RANCIF=4,ARANCIF=5)
  if (dep.model==0) stop("model must be COR, OR, RR, RANCIF, ARANCIF \n"); 
  if (dep.model<=4) rvdes <- matrix(0,1,1); 
###  if (dep.model==4 & !is.null(cif2) ) dep.model <- 6 ## two cause random cif
###  if (dep.model==6) stop("Different causes  under development \n"); 


  obj <- function(par) 
    { ## {{{

    if (is.null(par.func)==TRUE) {
	    Xtheta <- theta.des %*% par; DXtheta <- array(0,c(1,1,1));
    } else {
	  Xtheta <- c(); DXtheta <- array(0,c(length(times),antpers,dimpar));
	  if (is.null(dpar.func))
		  stop("Must also provide derivative of function wrt to parameters")
	  s <- 0
	  for (t in times) {
		  s <- 1+s
		  Xttheta <- par.func(par,t,theta.des)
		  if (length(Xttheta)!=antpers) stop("parfunc(par,t,theta.des) must length n"); 
		  Xtheta <- cbind(Xtheta,Xttheta)
		  Dttheta <- dpar.func(par,t,theta.des); 
		  if (dim(Dttheta)[1]!=antpers || dim(Dttheta)[2]!=dimpar) 
			  stop("dparfunc must return matrix n x dimpar when called on dparfunc(par,t,theta.des)"); 
		  DXtheta[s,,] <- Dttheta
	  }
###	  print(dim(Xtheta)); print(dim(DXtheta))
    }

      outl<-.Call("cor", ## {{{
                 itimes=times,iy=time,icause=cause,iCA1=cause1,iKMc=Gcx, 
                 iz=X,iest=matrix(est[,-1],ncol=ncol(est)-1),iZgamma=c(Zgamma),isemi=semi,izsem=Z, 
                 itheta=c(par),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),
		 ithetades=theta.des,
                 icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
                 iinverse=inverse,iCA2=cause2,
                 ix2=X2,isemi2=semi2,iest2=as.matrix(est2[,-1]),iZgamma2=c(Z2gamma2),
                 iflexfunc=flex.func,iiid=iid,isym=sym,iweights=weights, 
                 isamecens=as.numeric(same.cens),istabcens=as.numeric(stab.cens),
		 iKMtimes=Gctimes,isilent=silent,
                 icifmodel=cif.model,idepmodel=dep.model,
                 iestimator=estimator,ientryage=entry,icif1entry=cif1entry,
                 icif2entry=cif2entry,itrunkp=trunkp,irvdes=random.design
                 ### Biid,gamma.iid,time.pow, 
                 ### B2iid,gamma2.iid,body(htheta),body(dhtheta),new.env(),
		 ) 
      ## }}}

      attr(outl,"gradient") <-outl$score 
      if (oout==0) return(outl$ssf) else if (oout==1) return(sum(outl$score^2)) else return(outl)
    } ## }}}

  p <- theta
  iid <- 0;  ###no-iid representation for iterations
  if (score.method=="fisher.scoring") { ## {{{
    oout <- 2;  ### output control for obj
    if (Nit>0) 
    for (i in 1:Nit)
    {
        oout <- 2
        out <- obj(p)
	hess <- out$Dscore
	if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- hess 
        if (detail==1) {## {{{
          print(paste("Fisher-Scoring ===================: it=",i)); 
          cat("theta:");print(c(p))
          cat("score:");print(c(out$score)); 
	  cat("hess:"); print(hess); 
        }## }}}
        delta <- hessi %*% out$score 
	### for test purposes
###     oout <- 0;  ### output control for obj
###     score1 <- numDeriv::jacobian(obj,p)
###	hess1 <-  numDeriv::hessian(obj,p)
###	print(score1)
###	print(hess1)
	## do not update last iteration 
	if (i<Nit) p <- p-delta* step
	if (is.nan(sum(out$score))) break; 
        if (sum(abs(out$score))<0.00001) break; 
        if (max(delta)>20) { cat("NR increment > 20, lower step zize, increment= \n"); cat(delta); break; }
    }
    if (!is.nan(sum(p))) { ## {{{ iid decomposition
    oout <- 2
    theta <- p
    iid <- 1; 
    out <- obj(p) 
    score <- out$score
    hess <- out$Dscore
    } ## }}} 
    if (detail==1 & Nit==0) {## {{{
          print(paste("Fisher-Scoring ===================: final")); 
          cat("theta:");print(c(p))
          cat("score:");print(c(out$score)); 
	  cat("hess:"); print(hess); 
###	  oout <- 0; hess1 <-  numDeriv::hessian(obj,p); print(hess1)
    }## }}}
    if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- diag(nrow(hess))
###    score1 <- numDeriv::jacobian(obj,p)
    score1 <- score; 
    ## }}}
  } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlminb(theta,obj,control=control),error=function(x) NA)
    if (detail==1) print(opt); 
    iid <- 1; 
    hess <- numDeriv::hessian(obj,opt$par)
    score <- numDeriv::jacobian(obj,opt$par)
    hessi <- lava::Inverse(hess); 
    theta <- opt$par
    if (detail==1) cat("iid decomposition\n"); 
    oout <- 2; 
    out <- obj(opt$par)
    score1 <- out$score
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(obj,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1) cat("iid decomposition\n"); 
    oout <- 2; 
    out <- obj(opt$estimate)
    score1 <- out$score
  ## }}}
  }  else stop("score.methods = nlm nlminb fisher.scoring\n"); 

  theta.iid <- out$theta.iid %*% hessi
  if (is.null(par.func)) var.theta  <- t(theta.iid) %*% theta.iid else var.theta <- hessi
  var.theta  <- t(theta.iid) %*% theta.iid 
  if (is.null(par.func)) thetanames <- colnames(theta.des) else
  thetanames <- rep("R-func",dimpar)
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,
             theta.iid=theta.iid,score1=c(score1),thetanames=thetanames,
	     brierscore=out$ssf,p11=out$p11); 
  if (dep.model<=3) class(ud)<-"cor" 
  else if (dep.model==4) class(ud) <- "randomcif" 
  else if (dep.model==6) class(ud) <- "randomcif" 
  else if (dep.model==5) class(ud) <- "randomcifrv" 
  attr(ud, "Formula") <- formula
  attr(ud, "Clusters") <- clusters
  attr(ud,"cause1")<-cause1; attr(ud,"cause2")<-cause2
  attr(ud,"sym")<-sym; 
  attr(ud,"inverse")<-inverse; 
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  if (dep.model==4) attr(ud, "Type") <- "randomcif"
  if (dep.model==6) attr(ud, "Type") <- "randomcif"
  if (model=="COR") attr(ud, "Type") <- "cor"
  if (model=="RR") attr(ud, "Type") <- "RR"
  if (model=="OR") attr(ud, "Type") <- "OR-cif"
  return(ud);
} ## }}}

###mysolve <- function(A)
###{
###  ee <- eigen(A);
###  threshold <- 1e-12
###  idx <- ee$values>threshold
###  ee$values[idx] <- 1/ee$values[idx];
###  if (!all(idx))
###    ee$values[!idx] <- 0
###  V <- with(ee, vectors%*%diag(values)%*%t(vectors))
###  return(V)
###}

##' Fits a parametric model for the log-cross-odds-ratio for the 
##' predictive effect of for the cumulative incidence curves for \eqn{T_1} 
##' experiencing cause i given that \eqn{T_2} has experienced a cause k :
##' \deqn{
##' \log(COR(i|k))  =  h(\theta,z_1,i,z_2,k,t)=_{default}  \theta^T z =  
##' }
##' with the log cross odds ratio being 
##' \deqn{
##' COR(i|k) = 
##' \frac{O(T_1 \leq t,cause_1=i | T_2 \leq t,cause_2=k)}{
##' O(T_1 \leq t,cause_1=i)}  
##' }
##' the conditional odds divided by the unconditional odds, with the odds
##' being, respectively 
##' \deqn{
##' O(T_1 \leq t,cause_1=i | T_2 \leq t,cause_1=k) = 
##' \frac{
##' P_x(T_1 \leq t,cause_1=i | T_2 \leq t,cause_2=k)}{
##' P_x((T_1 \leq t,cause_1=i)^c | T_2 \leq t,cause_2=k)}
##' }
##' and 
##' \deqn{
##' O(T_1 \leq t,cause_1=i) = 
##' \frac{P_x(T_1 \leq t,cause_1=i )}{P_x((T_1 \leq t,cause_1=i)^c )}.
##' }
##' Here \eqn{B^c} is the complement event of \eqn{B},
##' \eqn{P_x} is the distribution given covariates 
##' (\eqn{x} are subject specific and \eqn{z} are cluster specific covariates), and 
##' \eqn{h()} is a function that is the simple identity 
##' \eqn{\theta^T z} by default.
##' 
##' The OR dependence measure is given by 
##' \deqn{
##' OR(i,k) = 
##' \log ( 
##' \frac{O(T_1 \leq t,cause_1=i | T_2 \leq t,cause_2=k)}{
##' O(T_1 \leq t,cause_1=i) | T_2 \leq t,cause_2=k)}  
##' }
##' This measure is numerically more stabile than the COR measure, and is symetric in i,k.
##' 
##' The RR dependence measure is given by 
##' \deqn{
##' RR(i,k) = 
##' \log ( 
##' \frac{P(T_1 \leq t,cause_1=i , T_2 \leq t,cause_2=k)}{
##' P(T_1 \leq t,cause_1=i) P(T_2 \leq t,cause_2=k)}  
##' }
##' This measure is numerically more stabile than the COR measure, and is symetric in i,k.
##' 
##' The model is fitted under symmetry (sym=1), i.e., such that it is assumed 
##' that \eqn{T_1} and \eqn{T_2} can be interchanged and leads to
##' the same cross-odd-ratio (i.e.
##' \eqn{COR(i|k) = COR(k|i))}, 
##' as would be expected for twins 
##' or without symmetry as might be the case with mothers and daughters (sym=0). 
##' 
##' \eqn{h()} may be specified as an R-function of the parameters, 
##' see example below, but the default is that it is simply \eqn{\theta^T z}.
##'
##' @title Cross-odds-ratio, OR or RR risk regression for competing risks
##' @aliases or.cif rr.cif
##' @param cif a model object from the comp.risk function with the 
##' marginal cumulative incidence of cause1, i.e., the event of interest, and whose
##' odds the comparision is compared to the conditional odds given cause2
##' @param data a data.frame with the variables.
##' @param cause specifies the causes  related to the death
##' times, the value cens.code is the censoring value. When missing it comes from marginal cif.
##' @param times time-vector that specifies the times used for the estimating euqations for the cross-odds-ratio estimation.
##' @param cause1 specificies the cause considered.
##' @param cause2 specificies the cause that is conditioned on.
##' @param cens.code specificies the code for the censoring if NULL then uses the one from the marginal cif model.
##' @param cens.model specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"
##' @param Nit number of iterations for Newton-Raphson algorithm.
##' @param detail if 0 no details are printed during iterations, if 1 details are given.
##' @param clusters specifies the cluster structure.
##' @param theta specifies starting values for the cross-odds-ratio parameters of the model.
##' @param theta.des specifies a regression design for the cross-odds-ratio parameters.
##' @param step specifies the step size for the Newton-Raphson algorithm.
##' @param sym specifies if symmetry is used in the model.
##' @param weights weights for estimating equations.
##' @param par.func parfunc
##' @param dpar.func dparfunc
##' @param dimpar dimpar
##' @param score.method "nlminb", can also use "fisher-scoring".
##' @param same.cens if true then censoring within clusters are assumed to be the same variable, default is independent censoring.
##' @param censoring.weights these probabilities are used for the bivariate censoring dist.
##' @param silent 1 to suppress output about convergence related issues.
##' @param ... Not used.
##' @return returns an object of type 'cor'. With the following arguments:
##' \item{theta}{estimate of proportional odds parameters of model.}
##' \item{var.theta}{variance for gamma.  }
##' \item{hess}{the derivative of the used score.}
##' \item{score}{scores at final stage.}
##' \item{score}{scores at final stage.}
##' \item{theta.iid}{matrix of iid decomposition of parametric effects.}
##' @author Thomas Scheike
##' @references
##' Cross odds ratio Modelling of dependence for
##' Multivariate Competing Risks Data, Scheike and Sun (2010), work in progress.
##' 
##' A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
##' Scheike, Zhang, Sun, Jensen (2010), Biometrika. 
##' @examples
##' data(multcif);
##' multcif$cause[multcif$cause==0] <- 2
##' zyg <- rep(rbinom(200,1,0.5),each=2)
##' theta.des <- model.matrix(~-1+factor(zyg))
##' 
##' times=seq(0.05,1,by=0.05) # to speed up computations use only these time-points
##' add<-comp.risk(Event(time,cause)~+1+cluster(id),data=multcif,cause=1,
##'                n.sim=0,times=times,model="fg",max.clust=NULL)
##' add2<-comp.risk(Event(time,cause)~+1+cluster(id),data=multcif,cause=2,
##'                n.sim=0,times=times,model="fg",max.clust=NULL)
##' 
##' out1<-cor.cif(add,data=multcif,cause1=1,cause2=1)
##' summary(out1)
##' 
##' out2<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta.des=theta.des)
##' summary(out2)
##' 
##' ##out3<-cor.cif(add,data=multcif,cause1=1,cause2=2,cif2=add2)
##' ##summary(out3)
##' ###########################################################
##' # investigating further models using parfunc and dparfunc
##' ###########################################################
##' \donttest{ ## Reduce Ex.Timings
##' set.seed(100)
##' prt<-simnordic.random(2000,cordz=2,cormz=5)
##' prt$status <-prt$cause
##' table(prt$status)
##' 
##' times <- seq(40,100,by=10)
##' cifmod <- comp.risk(Event(time,cause)~+1+cluster(id),data=prt,
##'                     cause=1,n.sim=0,
##'                     times=times,conservative=1,max.clust=NULL,model="fg")
##' theta.des <- model.matrix(~-1+factor(zyg),data=prt)
##' 
##' parfunc <- function(par,t,pardes)
##' {
##' par <- pardes %*% c(par[1],par[2]) +
##'        pardes %*% c( par[3]*(t-60)/12,par[4]*(t-60)/12)
##' par
##' }
##' head(parfunc(c(0.1,1,0.1,1),50,theta.des))
##' 
##' dparfunc <- function(par,t,pardes)
##' {
##' dpar <- cbind(pardes, t(t(pardes) * c( (t-60)/12,(t-60)/12)) )
##' dpar
##' }
##' head(dparfunc(c(0.1,1,0.1,1),50,theta.des))
##' 
##' names(prt)
##' or1 <- or.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,
##'               same.cens=TRUE,theta=c(0.6,1.1,0.1,0.1),
##'               par.func=parfunc,dpar.func=dparfunc,dimpar=4,
##'               score.method="fisher.scoring",detail=1)
##' summary(or1)
##' 
##'  cor1 <- cor.cif(cifmod,data=prt,cause1=1,cause2=1,theta.des=theta.des,
##'                  same.cens=TRUE,theta=c(0.5,1.0,0.1,0.1),
##'                  par.func=parfunc,dpar.func=dparfunc,dimpar=4,
##'                  control=list(trace=TRUE),detail=1)
##' summary(cor1)
##' 
##' ### piecewise contant OR model
##' gparfunc <- function(par,t,pardes)
##' {
##' 	cuts <- c(0,80,90,120)
##' 	grop <- diff(t<cuts)
##' paru  <- (pardes[,1]==1) * sum(grop*par[1:3]) +
##'     (pardes[,2]==1) * sum(grop*par[4:6])
##' paru
##' }
##' 
##' dgparfunc <- function(par,t,pardes)
##' {
##' 	cuts <- c(0,80,90,120)
##' 	grop <- diff(t<cuts)
##' par1 <- matrix(c(grop),nrow(pardes),length(grop),byrow=TRUE)
##' parmz <- par1* (pardes[,1]==1)
##' pardz <- (pardes[,2]==1) * par1
##' dpar <- cbind( parmz,pardz)
##' dpar
##' }
##' head(dgparfunc(rep(0.1,6),50,theta.des))
##' head(gparfunc(rep(0.1,6),50,theta.des))
##' 
##' or1g <- or.cif(cifmod,data=prt,cause1=1,cause2=1,
##'                theta.des=theta.des, same.cens=TRUE,
##'                par.func=gparfunc,dpar.func=dgparfunc,
##'                dimpar=6,score.method="fisher.scoring",detail=1)
##' summary(or1g)
##' names(or1g)
##' head(or1g$theta.iid)
##' }
##' @export
##' @keywords survival
cor.cif<-function(cif,data,cause=NULL,times=NULL,
                  cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
                  clusters=NULL, theta=NULL,theta.des=NULL,step=1,sym=0,weights=NULL, 
		  par.func=NULL,dpar.func=NULL,dimpar=NULL,
		  score.method="nlminb",same.cens=FALSE,censoring.weights=NULL,silent=1,...)
{ ## {{{
  fit <- dep.cif(cif=cif,data=data,cause=cause,model="COR",times=times,
                 cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
                 clusters=clusters,theta=theta,theta.des=theta.des,par.func=par.func,dpar.func=dpar.func,
		 dimpar=dimpar,
                 step=step,sym=sym,weights=weights,
                 score.method=score.method,same.cens=same.cens,censoring.weights=censoring.weights,silent=silent,...)
  fit$call <- match.call()
  fit
} ## }}}
##' @export
rr.cif<-function(cif,data,cause=NULL,cif2=NULL,times=NULL,
                 cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
                 clusters=NULL, theta=NULL,theta.des=NULL, step=1,sym=0,weights=NULL,
                 same.cens=FALSE,censoring.weights=NULL,silent=1,par.func=NULL,dpar.func=NULL,dimpar=NULL,
		 score.method="nlminb",entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,...)
{ ## {{{
  fit <- dep.cif(cif=cif,data=data,cause=cause,model="RR",cif2=cif2,times=times,
                 cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
                 clusters=clusters,theta=theta,theta.des=theta.des, step=step,sym=sym,weights=weights,
                 same.cens=same.cens,censoring.weights=censoring.weights,silent=silent,
		 par.func=par.func,dpar.func=dpar.func,dimpar=dimpar, score.method=score.method,
                 entry=entry,estimator=estimator,trunkp=trunkp,admin.cens=admin.cens,...)
  fit$call <- match.call()
  fit
} ## }}}
##' @export
or.cif<-function(cif,data,cause=NULL,cif2=NULL,times=NULL,
                 cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
                 clusters=NULL, theta=NULL,theta.des=NULL, step=1,sym=0, weights=NULL,
                 same.cens=FALSE,censoring.weights=NULL,silent=1,par.func=NULL,dpar.func=NULL,dimpar=NULL,
		 score.method="nlminb",entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,...)
{ ## {{{
  fit <- dep.cif(cif=cif,data=data,cause=cause,model="OR",cif2=cif2,times=times,
                 cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
                 clusters=clusters,theta=theta,theta.des=theta.des, step=step,sym=sym,weights=weights,
                 same.cens=same.cens,censoring.weights=censoring.weights,silent=silent,
		 par.func=par.func,dpar.func=dpar.func,dimpar=dimpar,
		 score.method=score.method,
                 entry=entry,estimator=estimator,trunkp=trunkp,admin.cens=admin.cens,...)
  fit$call <- match.call()
  fit
} ## }}}

##' Fits a random effects  model describing the dependence in the cumulative 
##' incidence curves for subjects within a cluster.  Given the gamma distributed
##' random effects it is assumed that the cumulative incidence curves are indpendent, and
##' that the marginal cumulative incidence curves are on the form
##' \deqn{
##' P(T \leq t, cause=1 | x,z) = P_1(t,x,z) = 1- exp( -x^T A(t) exp(z^T \beta))
##' }
##' We allow a regression structure for the random effects variances that may depend on
##' cluster covariates.
##'
##' @title Random effects model for competing risks data
##' @param cif a model object from the comp.risk function with the 
##' marginal cumulative incidence of cause2, i.e., the event that is conditioned on, and whose
##' odds the comparision is made with respect to
##' @param data a data.frame with the variables.
##' @param cause specifies the causes  related to the death
##' times, the value cens.code is the censoring value.
##' @param cif2 specificies model for cause2 if different from cause1.
##' @param cause1 cause of first coordinate.
##' @param cause2 cause of second coordinate.
##' @param cens.code specificies the code for the censoring if NULL then uses the one from the marginal cif model.
##' @param cens.model specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"
##' @param Nit number of iterations for Newton-Raphson algorithm.
##' @param detail if 0 no details are printed during iterations, if 1 details are given.
##' @param clusters specifies the cluster structure.
##' @param theta specifies starting values for the cross-odds-ratio parameters of the model.
##' @param theta.des specifies a regression design for the cross-odds-ratio parameters.
##' @param sym 1 for symmetry 0 otherwise
##' @param step specifies the step size for the Newton-Raphson algorith.m
##' @param same.cens if true then censoring within clusters are assumed to be the same variable, default is independent censoring.
##' @param exp.link if exp.link=1 then var is on log-scale.
##' @param score.method default uses "nlminb" optimzer, alternatively, use the "fisher-scoring" algorithm.
##' @param entry entry-age in case of delayed entry. Then two causes must be given.
##' @param trunkp gives probability of survival for delayed entry, and related to entry-ages given above.
##' @param ... extra arguments.
##' @return returns an object of type 'cor'. With the following arguments:
##' \item{theta}{estimate of proportional odds parameters of model.}
##' \item{var.theta}{variance for gamma.  }
##' \item{hess}{the derivative of the used score.}
##' \item{score}{scores at final stage.}
##' \item{score}{scores at final stage.}
##' \item{theta.iid}{matrix of iid decomposition of parametric effects.}
##' @export
##' @references
##' A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
##' Scheike, Zhang, Sun, Jensen (2010), Biometrika. 
##' 
##' Cross odds ratio Modelling of dependence for
##' Multivariate Competing Risks Data, Scheike and Sun (2012), work in progress.
##' @examples
##' \donttest{ ## Reduce Ex.Timings
##'  d <- simnordic.random(4000,delayed=TRUE,
##'        cordz=0.5,cormz=2,lam0=0.3,country=TRUE)
##'  times <- seq(50,90,by=10)
##'  add1<-comp.risk(Event(time,cause)~const(country)+cluster(id),data=d,
##'  times=times,cause=1,max.clust=NULL)
##' 
##'  ### making group indidcator 
##'  mm <- model.matrix(~-1+factor(zyg),d)
##' 
##'  out1<-random.cif(add1,data=d,cause1=1,cause2=1,theta=1,same.cens=TRUE)
##'  summary(out1)
##' 
##'  out2<-random.cif(add1,data=d,cause1=1,cause2=1,theta=1,
##' 		   theta.des=mm,same.cens=TRUE)
##'  summary(out2)
##' 
##' #########################################
##' ##### 2 different causes
##' #########################################
##' 
##'  add2<-comp.risk(Event(time,cause)~const(country)+cluster(id),data=d,
##'                   times=times,cause=2,max.clust=NULL)
##'  out3<-random.cif(add1,data=d,cause1=1,cause2=2,cif2=add2,sym=1,same.cens=TRUE)
##'  summary(out3) ## negative dependence
##' 
##'  out4<-random.cif(add1,data=d,cause1=1,cause2=2,cif2=add2,theta.des=mm,sym=1,same.cens=TRUE)
##'  summary(out4) ## negative dependence
##' }
##' @keywords survival
##' @author Thomas Scheike
random.cif<-function(cif,data,cause=NULL,cif2=NULL,
                     cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
                     clusters=NULL,theta=NULL,theta.des=NULL,sym=1,
                     step=1,same.cens=FALSE,exp.link=0,score.method="fisher.scoring",
                     entry=NULL,trunkp=1,...)
{ ## {{{
  fit <- dep.cif(cif,data=data,cause=cause,model="RANCIF",cif2=cif2,sym=sym,
     cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
     clusters=clusters,theta=theta,theta.des=theta.des,step=step,same.cens=same.cens,
     exp.link=exp.link,score.method=score.method,entry=entry,trunkp=trunkp,...)
  fit$call <- match.call()
  fit
} ## }}}

##'Fits a random effects  model describing the dependence in the cumulative 
##'incidence curves for subjects within a cluster.  Given the gamma distributed
##'random effects it is assumed that the cumulative incidence curves 
##'are indpendent, and that the marginal cumulative incidence curves 
##'are on additive form
##'\deqn{
##'P(T \leq t, cause=1 | x,z) = P_1(t,x,z) = 1- exp( -x^T A(t) - t z^T \beta)
##'}
##'
##'We allow a regression structure for the indenpendent gamma distributed 
##'random effects  and their variances that may depend on cluster covariates.
##'
##'random.design specificies the random effects for each subject within a cluster. This is
##'a matrix of 1's and 0's with dimension n x d.  With d random effects. 
##'For a cluster with two subjects, we let the random.design rows be 
##' \eqn{v_1} and \eqn{v_2}. 
##'Such that the random effects for subject 
##'1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
##'has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}. By construction
##'subjects 1's random effect are Gamma distributed with 
##'mean \eqn{\lambda_1/v_1^T \lambda}
##'and variance \eqn{\lambda_1/(v_1^T \lambda)^2}. Note that the random effect 
##'\eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
##'
##'The parameters \eqn{(\lambda_1,...,\lambda_d)}
##'are related to the parameters of the model
##'by a regression construction \eqn{pard} (d x k), that links the \eqn{d} 
##'\eqn{\lambda} parameters
##'with the (k) underlying \eqn{\theta} parameters 
##'\deqn{
##'\lambda = pard \theta 
##'}
##'
##' @title Additive Random effects model for competing risks data for polygenetic modelling
##'
##'
##' @param cif a model object from the comp.risk function with the 
##' marginal cumulative incidence of cause2, i.e., the event that is conditioned on, and whose
##' odds the comparision is made with respect to
##' @param data a data.frame with the variables.
##' @param cause specifies the causes  related to the death
##' times, the value cens.code is the censoring value.
##' @param cif2 specificies model for cause2 if different from cause1.
##' @param times time points
##' @param cause1 cause of first coordinate.
##' @param cause2 cause of second coordinate.
##' @param cens.code specificies the code for the censoring if NULL then uses the one from the marginal cif model.
##' @param cens.model specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"
##' @param Nit number of iterations for Newton-Raphson algorithm.
##' @param detail if 0 no details are printed during iterations, if 1 details are given.
##' @param clusters specifies the cluster structure.
##' @param theta specifies starting values for the cross-odds-ratio parameters of the model.
##' @param theta.des specifies a regression design for the cross-odds-ratio parameters.
##' @param step specifies the step size for the Newton-Raphson algorith.m
##' @param sym 1 for symmetri and 0 otherwise
##' @param weights weights for score equations.
##' @param same.cens if true then censoring within clusters are assumed to be the same variable, default is independent censoring.
##' @param censoring.weights Censoring probabilities
##' @param silent debug information 
##' @param exp.link if exp.link=1 then var is on log-scale.
##' @param score.method default uses "nlminb" optimzer, alternatively, use the "fisher-scoring" algorithm.
##' @param entry entry-age in case of delayed entry. Then two causes must be given.
##' @param estimator estimator
##' @param trunkp gives probability of survival for delayed entry, and related to entry-ages given above.
##' @param admin.cens Administrative censoring
##' @param random.design specifies a regression design of 0/1's for the random effects.
##' @param ... extra arguments.
##' @return returns an object of type 'random.cif'. With the following arguments:
##'\item{theta}{estimate of parameters of model.}
##'\item{var.theta}{variance for gamma.  }
##'\item{hess}{the derivative of the used score.}
##'\item{score}{scores at final stage.}
##'\item{theta.iid}{matrix of iid decomposition of parametric effects.}
##' @export
##' @references
##' A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
##' Scheike, Zhang, Sun, Jensen (2010), Biometrika.
##'
##' Cross odds ratio Modelling of dependence for
##' Multivariate Competing Risks Data, Scheike and Sun (2013), Biostatitistics.
##'
##' Scheike, Holst, Hjelmborg (2014),  LIDA,  
##' Estimating heritability for cause specific hazards based on twin data
##' @examples
##' \donttest{ ## Reduce Ex.Timings
##'  d <- simnordic.random(5000,delayed=TRUE,
##'        cordz=0.5,cormz=2,lam0=0.3,country=TRUE)
##'  times <- seq(50,90,by=10)
##'  addm<-comp.risk(Event(time,cause)~const(country)+cluster(id),data=d,
##'  times=times,cause=1,max.clust=NULL)
##' 
##'  ### making group indidcator 
##'  mm <- model.matrix(~-1+factor(zyg),d)
##' 
##'  out1m<-random.cif(addm,data=d,cause1=1,cause2=1,theta=1,
##' 		   theta.des=mm,same.cens=TRUE)
##'  summary(out1m)
##'  
##'  ## this model can also be formulated as a random effects model 
##'  ## but with different parameters
##'  out2m<-Grandom.cif(addm,data=d,cause1=1,cause2=1,
##' 		    theta=c(0.4,4),step=0.5,
##' 		    random.design=mm,same.cens=TRUE)
##'  summary(out2m)
##'  1/out2m$theta
##'  out1m$theta
##'  
##'  ####################################################################
##'  ################### ACE modelling of twin data #####################
##'  ####################################################################
##'  ### assume that zygbin gives the zygosity of mono and dizygotic twins
##'  ### 0 for mono and 1 for dizygotic twins. We now formulate and AC model
##'  zygbin <- d$zyg=="DZ"
##' 
##'  n <- nrow(d)
##'  ### random effects for each cluster
##'  des.rv <- cbind(mm,(zygbin==1)*rep(c(1,0)),(zygbin==1)*rep(c(0,1)),1)
##'  ### design making parameters half the variance for dizygotic components
##'  pardes <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0), c(0,1))
##' 
##'  outacem <-Grandom.cif(addm,data=d,cause1=1,cause2=1,
##' 		same.cens=TRUE,theta=c(0.7,-0.3),
##'             step=1.0,theta.des=pardes,random.design=des.rv)
##'  summary(outacem)
##'  ## genetic variance is 
##'  exp(outacem$theta[1])/sum(exp(outacem$theta))^2
##' }
##' @keywords survival
##' @author Thomas Scheike
##' @export
Grandom.cif<-function(cif,data,cause=NULL,cif2=NULL,times=NULL,
cause1=1,cause2=1,cens.code=NULL,cens.model="KM",Nit=40,detail=0,
clusters=NULL, theta=NULL,theta.des=NULL, weights=NULL, step=1,sym=0,
same.cens=FALSE,censoring.weights=NULL,silent=1,exp.link=0,score.method="fisher.scoring",
entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,random.design=NULL,...)
{ ## {{{
fit <- dep.cif(cif,data=data,cause=cause,model="ARANCIF",cif2=cif2,times=times,
         cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
         clusters=clusters,theta=theta,theta.des=theta.des,step=step,sym=sym,weights=weights,
         same.cens=same.cens,censoring.weights=censoring.weights,silent=silent,exp.link=exp.link,
	 score.method=score.method,entry=entry,estimator=estimator,
	 random.design=random.design,trunkp=trunkp,admin.cens=admin.cens,...)
    fit$call <- match.call()
    fit
} ## }}}

##' @export
print.summary.cor <- function(x,digits=3,...)
{ ## {{{
  if (x$type=="cor") {
    cat("Cross odds ratio dependence for competing risks\n\n")
    cat("Effect of cause1=",x$cause1," on cause2=",x$cause2,
        " under symmetry=",x$sym,fill=TRUE,sep="")
  } else if (x$type=="RR") {
    cat("Ratio of joint and product of marginals for competing risks\n\n")
    cat("Ratio of cumulative incidence for cause1=",x$cause1," and cause2=",x$cause2,sep=" ")
  } else if (x$type=="OR-cif") {
    cat("OR for dependence for competing risks\n\n")
    cat("OR of cumulative incidence for cause1=",x$cause1," and cause2=",x$cause2,sep=" ")
  }
  cat("\n")
  prmatrix(signif(x$estimates,digits))
  cat("\n")
  if (!is.null(x$marg)) {
    cat(paste("Marginal cumulative incidencen",signif(x$marg,digits),"\n"))
    prmatrix(signif(x$casewise,digits))
    prmatrix(signif(x$concordance,digits))
    cat("\n")
  }
  invisible(x)
} ## }}}

##' Computes concordance and casewise concordance for dependence models for competing risks 
##' models of the type cor.cif, rr.cif or or.cif for the given cumulative incidences and the different dependence
##' measures in the object.
##'
##' @title Summary for dependence models for competing risks
##' @param object object from cor.cif rr.cif or or.cif for dependence between competing risks data for two causes.
##' @param marg.cif a number that gives the cumulative incidence in one time point for which concordance and 
##' casewise concordance are computed.
##' @param marg.cif2 the cumulative incidence for cause 2 for concordance and 
##' casewise concordance are computed. Default is that it is the same as marg.cif.
##' @param digits digits in output.
##' @param ... Additional arguments.
##' @return prints summary for dependence model. 
##' \item{casewise}{gives casewise concordance that is, probability of cause 2 (related to cif2) given that cause 1 (related to cif1)
##' 	has occured.}
##' \item{concordance}{gives concordance that is, probability of cause 2 (related to cif2) and cause 1 (related to cif1).}
##' \item{cif1}{cumulative incidence for cause1.}
##' \item{cif2}{cumulative incidence for cause1.}
##' @references
##' Cross odds ratio Modelling of dependence for
##' Multivariate Competing Risks Data, Scheike and Sun (2012), Biostatistics to appear.
##' 
##' A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
##' Scheike, Zhang, Sun, Jensen (2010), Biometrika. 
##' @author Thomas Scheike
##' @keywords survival
##' @examples
##' data(multcif) # simulated data 
##' multcif$cause[multcif$cause==0] <- 2
##' 
##' times=seq(0.1,3,by=0.1) # to speed up computations use only these time-points
##' add<-comp.risk(Event(time,cause)~const(X)+cluster(id),data=multcif,
##'                n.sim=0,times=times,cause=1)
##' ###
##' out1<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta=log(2+1))
##' summary(out1)
##' 
##' pad <- predict(add,X=1,Z=0,se=0,uniform=0)
##' summary(out1,marg.cif=pad)
##' @method summary cor
##' @export
summary.cor <- function(object,marg.cif=NULL,marg.cif2=NULL,digits=3,...) { ## {{{
  if (!(attr(object,"class") %in% c("cor","randomcif"))) stop("Must be a cor.cif or randomcif object")
  if (sum(abs(object$score))>0.001) warning("WARNING: check score for convergence\n")

  coefs <- coef(object,...) 

  ocasewise <- oconcordance <- NULL
  if (is.null(marg.cif)==FALSE) { ## {{{ 
    time <- marg.cif$time
    marg.cif <- marg.cif$P1
    if (is.null(marg.cif2)==FALSE) marg.cif2 <- marg.cif2$P1 else {
    if (attr(object,"cause2")==attr(object,"cause1")) marg.cif2 <- marg.cif 
    else stop("causes not the same and second marginal cif not given\n"); 
    }
    pmarg.cif <- marg.cif*marg.cif2

    thetav <- coefs[,1]
    thetavl <- thetav-1.96*coefs[,2]
    thetavu <- thetav+1.96*coefs[,2]
    namev <- object$thetanames
    if (is.null(namev)) namev <- paste("name",1:length(thetav),sep="")

    ocasewise <- oconcordance <- list()
    k <- 0
    for (theta in thetav) { 
       k <- k+1
       thetal <- thetavl[k]; 
       thetau <- thetavu[k]
    if (attr(object,"Type")=="cor") { ## {{{
      concordance <- exp(theta)*pmarg.cif/((1-marg.cif)+exp(theta)*marg.cif)
      conclower  <- exp(thetal)*pmarg.cif/((1-marg.cif)+exp(thetal)*marg.cif)
      concup  <- exp(thetau)*pmarg.cif/((1-marg.cif)+exp(thetau)*marg.cif)
      casewise <- concordance/marg.cif
      caselower  <- conclower/marg.cif
      caseup     <- concup/marg.cif
    } else if (attr(object,"Type")=="RR") {
      casewise<- exp(coefs[,1])*c(marg.cif)
      concordance <- exp(coefs[,1])*pmarg.cif
      caselower  <- marg.cif*exp(thetal)
      caseup     <- marg.cif*exp(thetau)
      conclower  <- pmarg.cif* exp(thetal)
      concup     <- pmarg.cif*exp(thetau)
    } else if (attr(object,"Type")=="OR-cif") {
      casewise<- plack.cif2(marg.cif,marg.cif,c(theta))/marg.cif
      concordance <- plack.cif2(marg.cif,marg.cif,c(theta))
      caselower  <- plack.cif2(marg.cif,marg.cif,thetal)/marg.cif
      caseup  <- plack.cif2(marg.cif,marg.cif,thetau)/marg.cif
      conclower  <- plack.cif2(marg.cif,marg.cif,thetal)
      concup     <- plack.cif2(marg.cif,marg.cif,thetau)
    } else if (attr(object,"Type")=="randomcif") {
	    theta <- 1/theta
	    thetal <- 1/thetal
	    thetau <- 1/thetau
      lam <- 1-marg.cif
      p11<- 1-lam -lam +lap(theta,2*ilap(theta, lam))
      concordance <- p11 
      conclower  <- 1-lam -lam +lap(thetal,2*ilap(thetal, lam))
      concup     <- 1-lam -lam +lap(thetau,2*ilap(thetau, lam))
      casewise <- concordance/marg.cif
      caselower  <- conclower/marg.cif
      caseup     <- concup/marg.cif
    } ## }}} 
    outcase <- cbind(time,casewise,caselower,caseup)
    outconc <- cbind(time,concordance,conclower,concup)
###    rownames(outcase) <- rownames(outconc)  <-  rownames(coefs)
    colnames(outcase) <- c("time","casewise concordance","2.5 %","97.5%")
    colnames(outconc) <- c("time","concordance","2.5 %","97.5%")
	  if (length(thetav)==1) {
	     ocasewise <- outcase
	     oconcordance <- outconc
	  } else {
	     ocasewise[[k]] <- outcase
	     oconcordance[[k]] <- outconc
	  }
  } 

  if (length(thetav)>1) {
  if (length(ocasewise)==length(namev)) {
    names(ocasewise) <- namev
    names(oconcordance) <- namev
  }
  }

  } ## }}} 

  res <- list(casewise=ocasewise,concordance=oconcordance,estimates=coefs,
	      marg.cif=marg.cif, marg.cif2=marg.cif2,type=attr(object,"Type"),
	      sym=attr(object,"sym"),cause1=attr(object,"cause1"),cause2=attr(object,"cause2"))
  class(res) <- "summary.cor"
  res
} ## }}}

##' @export
coef.cor <- function(object,...)
{ ## {{{
  res <- cbind(object$theta, diag(object$var.theta)^0.5)
  se<-diag(object$var.theta)^0.5
  wald <- object$theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  cor<-exp(object$theta)
  res <- as.matrix(cbind(res, wald, waldp,cor,se*cor))
  if (attr(object,"Type")=="cor") 
    colnames(res) <- c("log-Coef.", "SE", "z", "P-val","Cross odds ratio","SE")
  else colnames(res) <- c("log-ratio Coef.", "SE", "z", "P-val","Ratio","SE")
  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames

  return(res)
} ## }}}

##' @export
print.cor<-function(x,digits=3,...)
{ ## {{{
  print(x$call); 
  cat("\n")
  print(summary(x)); 
} ## }}}

##' Concordance
##'
##' @title Concordance Computes concordance and casewise concordance
##' @param object Output from the cor.cif, rr.cif or or.cif function
##' @param cif1 Marginal cumulative incidence
##' @param cif2 Marginal cumulative incidence of other cause (cause2) if  it is different from cause1
##' @param messages  To print messages
##' @param model  Specfifies wich model that is considered if object not given.
##' @param coefs Specfifies  dependence parameters if object is not given.
##' @param ...  Extra arguments, not used.
##' @author Thomas Scheike
##' @export
concordance <- function(object,cif1,cif2=NULL,messages=TRUE,model=NULL,coefs=NULL,...)
{ ## {{{

  if (is.null(model)) { 
    if (!inherits(object, "cor")) stop("Must be a rr.cif, cor.cif or or.cif object")
    model <- attr(object,"Type")
  } 
  if (is.null(coefs)) coefs <- coef(object)
  if (is.null(coefs)) stop("Must give dependence parameters\n"); 

  if (!is.null(object)) { 
    cause1 <-  attr(object,"cause1");
    cause2 <-  attr(object,"cause2"); 
  } else cause1 <- cause2 <- 1

  if (is.null(cif2)==TRUE) cif2 <- cif1; 

  if (messages) {
    if (model=="cor") { ## {{{
      message("Cross odds ratio dependence for competing risks\n\n")
      message("Odds of cause1=",cause1," given cause2=",cause2," relative to Odds of cause1=",cause1,"\n",fill=TRUE,sep="")
    } else if (model=="RR") {
      message("Ratio of joint and product of marginals for competing risks\n\n")
      message("Ratio of cumulative incidence for cause1=",cause1," and cause2=",cause2,sep=" ")
    } else if (model=="OR-cif") {
      message("OR for dependence for competing risks\n\n")
      message("OR of cumulative incidence for cause1=",cause1," and cause2=",cause2,sep=" ")
    } ## }}}
  }

  out <- list()
  for (k in 1:nrow(coefs)) {
    ## {{{
    if (model=="cor") {
      concordance <- exp(coefs[k,1])*cif1*cif2/((1-cif1)+exp(coefs[k,1])*cif1)
      conclower  <- exp(coefs[k,1]-1.96*coefs[k,2])*cif1*cif2/((1-cif1)+exp(coefs[k,1]-1.96*coefs[k,2])*cif1)
      concup  <- exp(coefs[k,1]+1.96*coefs[k,2])*cif1*cif2/((1-cif1)+exp(coefs[k,1]+1.96*coefs[k,2])*cif1)
      casewise <- concordance/cif1
      caselower  <- conclower/cif1
      caseup     <- concup/cif1
    } else if (model=="RR") {
      casewise<- exp(coefs[k,1])*c(cif2)
      concordance <- exp(coefs[k,1])*cif1*cif2
      caselower  <- cif2*exp(coefs[k,1]-1.96*coefs[k,2])
      caseup     <- cif2*exp(coefs[k,1]+1.96*coefs[k,2])
      conclower  <- cif1*cif2* exp(coefs[k,1]-1.96*coefs[k,2])
      concup     <- cif1*cif2*exp(coefs[k,1]+1.96*coefs[k,2])
    } else if (model=="OR-cif") {
      thetal <-  coefs[k,1]-1.96*coefs[k,2]
      thetau <-  coefs[k,1]+1.96*coefs[k,2]
      casewise<- plack.cif2(cif1,cif2,c(coefs[k,1]))/cif1
      concordance <- plack.cif2(cif1,cif2,c(coefs[k,1]))
      caselower  <- plack.cif2(cif1,cif2,thetal)/cif1
      caseup  <- plack.cif2(cif1,cif2,thetau)/cif1
      conclower  <- plack.cif2(cif1,cif2,thetal)
      concup     <- plack.cif2(cif1,cif2,thetau)
    }

    outcase <- cbind(c(casewise),c(caselower),c(caseup))
    outconc <- cbind(c(concordance),c(conclower),c(concup))
    colnames(outcase) <- c("casewise concordance","2.5 %","97.5%")
    colnames(outconc) <- c("concordance","2.5 %","97.5%")
    ## }}}

    out[[k]] <- list(concordance=outconc,casewise=outcase)
    names(out)[k] <- rownames(coefs)[k]
###k <- k+1
  }

  return(out)
} ## }}}

##' .. content for description (no empty lines) ..
##'
##' @title plack Computes concordance for or.cif based model, that is Plackett random effects model 
##' @aliases plack.cif2
##' @export
##' @param cif1 Cumulative incidence of first argument.
##' @param cif2 Cumulative incidence of second argument.
##' @param object or.cif object with dependence parameters.
##' @author Thomas Scheike
plack.cif <- function(cif1,cif2,object) 
{ ## {{{
  coefs <- coef(object)
  theta <- exp(object$theta); 
  cif1 <- c(cif1); cif2 <- c(cif2)
  cifs=cif1+cif2; 

  valn=2*(theta-1); 
  val1=(1+(theta-1)*(cifs))-( ((1+(theta-1)*cifs))^2-4*cif1*cif2*theta*(theta-1))^0.5; 
  vali=cif1*cif2;
  valr <- vali;
  valr[valn!=0] <- val1/valn; 

  valr <- matrix(valr,length(c(theta)),1)
  rownames(valr)=colnames(coefs)
  return(valr); 
} ## }}}

##' @export
plack.cif2 <- function(cif1,cif2,theta) 
{ ## {{{
  theta <- exp(c(theta))
  cif1 <- c(cif1); cif2 <- c(cif2)
  cifs=cif1+cif2; 

  valn=2*(theta-1); 
  val1=(1+(theta-1)*(cifs))-( ((1+(theta-1)*cifs))^2-4*cif1*cif2*theta*(theta-1))^0.5; 
  vali=cif1*cif2;

  valr <- vali;
  valr[valn!=0] <- val1/valn; 

  return(valr); 
} ## }}}

##' @export
summary.randomcif<-function (object, ...) 
{ ## {{{
  if (!inherits(object, "randomcif")) 
    stop("Must be a random.cif  object")
  cat("Random effect variance for variation due to clusters\n\n")
  cat("Cause", attr(object, "cause1"), "and cause", attr(object, 
                                                         "cause2"), fill = TRUE)
  cat("\n")
  if (sum(abs(object$score)) > 0.001) cat("WARNING: check score for convergence")
  cat("\n")
  coef.randomcif(object, ...)
} ## }}}

##' @export
coef.randomcif<- function (object, digits = 3, ...) 
{ ## {{{
  res <- cbind(object$theta, diag(object$var.theta)^0.5)
  se <- diag(object$var.theta)^0.5
  wald <- object$theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  cor <- object$theta + 1
  res <- as.matrix(cbind(res, wald, waldp, cor, se))
  colnames(res) <- c("Coef.", "SE", "z", "P-val", "Cross odds ratio", "SE")
  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames

###  prmatrix(signif(res, digits))
  return(res)
} ## }}}

##' @export
print.randomcif<- function (x , digits = 3, ...) 
{ ## {{{
} ## }}}

##' @export
summary.randomcifrv<-function (object, ...) 
{ ## {{{
    if (!inherits(object, "randomcifrv")) 
        stop("Must be a random.cifrv  object")
    cat("Random effect parameters for additive gamma random effects \n\n")
    cat("Cause", attr(object, "cause1"), "and cause", attr(object, 
        "cause2"), fill = TRUE)
    cat("\n")
    if (sum(abs(object$score)) > 1e-06) 
        cat("WARNING: check score for convergence")
    cat("\n")
    coef.randomcifrv(object, ...)
} ## }}}

##' @export
coef.randomcifrv<- function (object, digits = 3, ...) 
{ ## {{{
    if (attr(object,"inverse")==1) elog <- 1 else elog  <- 0; 
    if (elog==1) theta <- exp(object$theta) else theta <- object$theta 
    se <- diag(object$var.theta)^0.5
    res <- cbind(object$theta, se)
    wald <- object$theta/se
    waldp <- (1 - pnorm(abs(wald))) * 2
    res <- as.matrix(cbind(res, wald, waldp))
    if (elog==0)  colnames(res) <- c("Coef.", "SE", "z", "P-val")
    if (elog==1) res <- cbind(res,exp(object$theta), exp(object$theta)^2*se)  
    if (elog==1) colnames(res) <- c("log-parameter","SE","z","P-val","exp(theta)","SE")

   if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
   prmatrix(signif(res, digits))

    cat("\n\n Random effect variances for gamma random effects \n\n")
    varpar <- theta/sum(theta)^2 
    res <- as.matrix(varpar); 
    if (elog==0)  { var.theta <-   object$var.theta; 
                    df <- 0*var.theta; 
                    for (i in 1:nrow(var.theta))
                    df[i,] <- -theta[i]*2*theta; 
		    diag(df) <- diag(df)+sum(theta)^2
		    df <- df/sum(theta)^4
		    var.varpar <- df %*% var.theta %*% df
                  }
    if (elog==1)  { 
	            var.theta <-   object$var.theta; 
                    var.varpar <- var.theta
                  }
    res <- cbind(res,diag(var.varpar)^.5)  
    colnames(res) <- c("variance","SE")
    if (is.null((rownames(res))) == TRUE) rownames(res) <- rep(" ", nrow(res))
    prmatrix(signif(res, digits))

} ## }}}

##' @export
print.randomcifrv<- function (x , digits = 3, ...) 
{ ## {{{
 summary(x, ...)
} ## }}}

## summary.cor<-function(object,digits=3,marg.cif=NULL,...)
## { ## {{{
##   if (!inherits(object, "cor")) stop("Must be a cor.cif  object")
##   if (sum(abs(object$score))>0.001) warning("WARNING: check score for convergence\n")
##   coefs <- coef.cor(object,...);

##   outcase <- outconc <- NULL
##   if (is.null(marg.cif)==FALSE) {
##     marg.cif <- max(marg.cif)
##     ## {{{
##     if (attr(object,"Type")=="cor") {
##       concordance <- exp(coefs[,1])*marg.cif^2/((1-marg.cif)+exp(coefs[,1])*marg.cif)
##       conclower  <- exp(coefs[,1]-1.96*coefs[,2])*marg.cif^2/((1-marg.cif)+exp(coefs[,1]-1.96*coefs[,2])*marg.cif)
##       concup  <- exp(coefs[,1]+1.96*coefs[,2])*marg.cif^2/((1-marg.cif)+exp(coefs[,1]+1.96*coefs[,2])*marg.cif)
##       casewise <- concordance/marg.cif
##       caselower  <- conclower/marg.cif
##       caseup     <- concup/marg.cif
##     } else if (attr(object,"Type")=="RR") {
##       casewise<- exp(coefs[,1])*c(marg.cif)
##       concordance <- exp(coefs[,1])*marg.cif^2
##       caselower  <- marg.cif*exp(coefs[,1]-1.96*coefs[,2])
##       caseup     <- marg.cif*exp(coefs[,1]+1.96*coefs[,2])
##       conclower  <- marg.cif^2* exp(coefs[,1]-1.96*coefs[,2])
##       concup     <- marg.cif^2*exp(coefs[,1]+1.96*coefs[,2])
##     } else if (attr(object,"Type")=="OR-cif") {
##       thetal <-  coefs[,1]-1.96*coefs[,2]
##       thetau <-  coefs[,1]+1.96*coefs[,2]
##       casewise<- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))/marg.cif
##       concordance <- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))
##       caselower  <- plack.cif2(marg.cif,marg.cif,thetal)/marg.cif
##       caseup  <- plack.cif2(marg.cif,marg.cif,thetau)/marg.cif
##       conclower  <- plack.cif2(marg.cif,marg.cif,thetal)
##       concup     <- plack.cif2(marg.cif,marg.cif,thetau)
##     }
##     outcase <- cbind(casewise,caselower,caseup)
##     outconc <- cbind(concordance,conclower,concup)
##     rownames(outcase) <- rownames(outconc)  <-  rownames(coefs)
##     colnames(outcase) <- c("casewise concordance","2.5 %","97.5%")
##     colnames(outconc) <- c("concordance","2.5 %","97.5%")
##   }
##   ## }}}
##   res <- list(casewise=outcase,concordance=outconc,estimates=coefs,marg=marg.cif,type=attr(object,"Type"),sym=attr(object,"sym"),cause1=attr(object,"cause1"),cause2=attr(object,"cause2"))
##   class(res) <- "summary.cor"
##   res
## } ## }}}

## coef.cor<-function(object,...)
## { ## {{{
##   res <- cbind(object$theta, diag(object$var.theta)^0.5)
##   se<-diag(object$var.theta)^0.5
##   wald <- object$theta/se
##   waldp <- (1 - pnorm(abs(wald))) * 2
##   cor<-exp(object$theta)
##   res <- as.matrix(cbind(res, wald, waldp,cor,se*cor))
##   if (attr(object,"Type")=="cor") 
##     colnames(res) <- c("log-Coef.", "SE", "z", "P-val","Cross odds ratio","SE")
##   else colnames(res) <- c("log-ratio Coef.", "SE", "z", "P-val","Ratio","SE")
##   if (is.null((rownames(res)))==TRUE) rownames(res)<-rep(" ",nrow(res))

##   return(res)
## } ## }}}
