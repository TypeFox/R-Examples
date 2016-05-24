comp.risk<-function(formula,data=sys.parent(),cause,times=NULL,Nit=50,clusters=NULL,est=NULL,
		    fix.gamma=0,gamma=0,n.sim=0,weighted=0,model="fg",detail=0,interval=0.01,resample.iid=1,
                    cens.model="KM",cens.formula=NULL,time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6,
                    weights=NULL,max.clust=1000,n.times=50,first.time.p=0.05,estimator=1,
		    trunc.p=NULL,cens.weights=NULL,admin.cens=NULL,conservative=1,monotone=0,step=NULL) 
    # {{{
{
    if (!missing(cause)){
       if (length(cause)!=1) stop("Argument cause has new meaning since 
   timereg version 1.8.4., it now specifies the cause of interest, see help(comp.risk) for details.")
   } 

    ## {{{
    # trans=1 P_1=1-exp( - ( x' b(b)+ z' gam t) ), 
    # trans=2 P_1=1-exp(-exp(x a(t)+ z` b )  Fine-Gray model, with baseline exp(x a(t)) 
    # trans=3 P_1= exp(x a(t)+ z` b)/( exp(x a(t) + z' b) +1 );  logistic
    # trans=4 P_1=exp( ( x' b(b)+ z' gam ) ), 
    # trans=5 P_1= (x' b(t)) exp( z' gam ), 
    # trans=6 P_1=1-exp(-(x a(t)) exp(z` b )) Fine-Gray model, with baseline x a(t) 
    # trans=7 P_1= (x a(t)) exp( z` b)/( (x a(t) ) exp(z' b) +1 ); logistic2
    trans <- switch(model,additive=1,prop=2,logistic=3,rcif=4,rcif2=5,fg=6,logistic2=7)
    ###  if (model=="additive")  trans<-1; if (model=="prop")      trans<-2; if (model=="logistic")  trans<-3; 
    ###  if (model=="rcif")      trans<-4; if (model=="rcif2")     trans<-5; if (model=="fg")        trans<-6; 
    ###  if (model=="logistic2") trans<-7; 
    line <- 0
    cause.call <- causeS <- cause
    m<-match.call(expand.dots=FALSE);
    m$gamma<-m$times<-m$n.times<-m$cause<-m$Nit<-m$weighted<-m$n.sim<-
             m$model<-m$detail<- m$cens.model<-m$time.pow<-m$silent<- m$step <- 
             m$cens.formula <- m$interval<- m$clusters<-m$resample.iid<- m$monotone <- 
             m$time.pow.test<-m$conv<- m$weights  <- m$max.clust <- m$first.time.p<- m$trunc.p <- 
             m$cens.weights <- m$admin.cens <- m$fix.gamma <- m$est  <- m$conservative <- m$estimator <- NULL
  
    if ((trans==2 || trans==3 || trans==7) && is.null(step)) step <- 0.5
    if (is.null(step)) step <- 1

    special <- c("const","cluster")
    if (missing(data)) {
        Terms <- terms(formula, special)
    }  else {
        Terms <- terms(formula, special, data = data)
    }
    m$formula <- Terms

    if (substr(as.character(m$formula)[2],1,4)=="Hist") {
       stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as 
       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
    }

    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (NROW(m) == 0) stop("No (non-missing) observations")
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    event.history <- model.extract(m, "response")

  if (class(event.history)!="Event"){
       stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as 
       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
  }

   model.type <- "competing.risks"

   ## {{{ Event stuff
    cens.code <- attr(event.history,"cens.code")
    if (ncol(event.history)==2) {
	time2 <- eventtime <- event.history[,1]
	status <- delta  <- event.history[,2]
	entrytime <- rep(0,length(time2))
	left <- 0
    } else {
	time2 <- eventtime <- event.history[,2]
	status <- delta  <- event.history[,3]
	entrytime <- event.history[,1]
	left <- 1
	if (max(entrytime)==0) left <- 0
    }
    event <- (status==cause)
    if (sum(event)==0) stop("No events of interest in data\n"); 

    ## }}} 

  if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;

  if (nrow(X)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 

  if (is.diag(t(X) %*% X)==TRUE) stratum <- 1 else stratum <- 0; 

  ## {{{ cluster set up
  if(is.null(clusters)){ clusters <- des$clusters}
  if(is.null(clusters)){
    cluster.call<-clusters; 
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    cluster.call<-clusters; 
    antclust <- length(unique(clusters))
    clusters <- as.integer(factor(clusters,labels=1:antclust))-1
  }

    coarse.clust <- FALSE; 
    if ((!is.null(max.clust))) if (max.clust< antclust) {
        coarse.clust <- TRUE
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
	antclust <- max.clust    
    }                                                         
    ## }}} 

    pxz <-px+pz;
    
    if (is.null(times)) {
        timesc<-sort(unique(eventtime[event==1])); 
        if (!is.null(n.times)) {
            if (length(timesc)> n.times) times <- quantile(timesc,prob=seq(first.time.p,1,length=n.times)) 
            else times <- timesc
        } else {times<-timesc; times<-times[times> quantile(timesc,prob=first.time.p)]; }
    } else times <- sort(times); 

  n<-nrow(X); ntimes<-length(times);
  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
  if (is.null(weights)==TRUE) weights <- rep(1,n); 
  ## }}}

  ## {{{ censoring and estimator 
  if (!is.null(admin.cens)) estimator  <- 3;
  Gcxe <- 1;  
  ordertime <- order(eventtime); 
  ###dcumhazcens <- rep(0,n); 

    if (estimator==1 || estimator==2) {
        if (is.null(cens.weights)) { ## {{{ censoring model stuff with possible truncation
            if (cens.model=="KM") { ## {{{
	        if (left==1) ud.cens<-survfit(Surv(entrytime,eventtime,delta==cens.code)~+1) else 
		ud.cens<-survfit(Surv(eventtime,delta==cens.code)~+1)
                Gfit<-cbind(ud.cens$time,ud.cens$surv)
                Gfit<-rbind(c(0,1),Gfit); 
                Gcx<-Cpred(Gfit,eventtime,strict=TRUE)[,2];
                Gcxe<-Cpred(Gfit,entrytime,strict=TRUE)[,2];
		### strictly before, but starts in 1. 
		Gcxe[Gcxe==0] <- 1
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; ## }}}
            } else if (cens.model=="stratKM") { ## {{{
	        XZ <- model.matrix(cens.formula,data=data); 
	        strata <- as.factor(XZ)
		Gcx <- pred.stratKM(data,time=eventtime,cause=delta,strata=strata)
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gctimes<-Cpred(Gfit,times)[,2]; ## }}}
            } else if (cens.model=="cox") { ## {{{
                if (!is.null(cens.formula)) { 
		      XZ <- model.matrix(cens.formula,data=data); 
                      if (sum(XZ[,1])==nrow(XZ)) XZ <- as.matrix(XZ[,-1])
                } else {
                       if (npar==TRUE) XZ<-X[,-1] else XZ <-cbind(X,Z)[,-1];
                }
		if (left==1) ud.cens<-coxph(Surv(entrytime,eventtime,delta==cens.code)~XZ)                
		else ud.cens<-coxph(Surv(eventtime,delta==cens.code)~XZ)                
		baseout <- basehaz(ud.cens,centered=FALSE); 
		baseout <- cbind(baseout$time,baseout$hazard)
		Gcx<-Cpred(baseout,eventtime,strict=TRUE)[,2];
		Gcxe<-Cpred(baseout,entrytime,strict=TRUE)[,2];
		Gcxe[Gcxe==0] <- 1
		RR<-exp(as.matrix(XZ) %*% coef(ud.cens))
		Gcx<-exp(-Gcx*RR)
		Gcxe<-exp(-Gcxe*RR)
		Gfit<-rbind(c(0,1),cbind(eventtime,Gcx)); 
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
		Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; 
                ## }}}
            } else if (cens.model=="aalen") {  ## {{{
                if (!is.null(cens.formula)) { 
			XZ <- model.matrix(cens.formula,data=data); 
                 } else {
                      if (npar==TRUE) XZ <-X else XZ <-cbind(X,Z);
                }
	        if (left==1) ud.cens<-aalen(Surv(entrytime,eventtime,delta==cens.code)~-1+XZ+cluster(clusters),
			       n.sim=0,residuals=0,robust=0,silent=1)
	        else ud.cens<-aalen(Surv(eventtime,delta==cens.code)~-1+XZ+cluster(clusters),
			       n.sim=0,residuals=0,robust=0,silent=1); 
                Gcx <- Cpred(ud.cens$cum,eventtime,strict=TRUE)[,-1];
                Gcx<-exp(-apply(Gcx*XZ,1,sum))
                Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-1
                Gcxe <- Cpred(ud.cens$cum,entrytime,strict=TRUE)[,2];
		Gcxe[Gcxe==0] <- 1
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gfit<-rbind(c(0,1),cbind(eventtime,Gcx)); 
                Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; ## }}}
            } else  stop('Unknown censoring model') 
            cens.weights <- Gcx
            if ((min(Gcx[event==1])< 0.00001) && (silent==0)) { 
                cat("Censoring dist. approx zero for some points, summary cens:\n");
                print(summary(Gcx)) 
            }
            ## }}}
        } else { 
            if (length(cens.weights)!=n) stop("censoring weights must have length equal to nrow in data\n");  
            Gcx <- cens.weights
	    ### for left truncation specification
            ord2 <- order(time2)
            Gctimes <- Cpred(cbind(time2[ord2],weights[ord2]),times)
        }
    } else { ## estimator==3 admin.cens 
        if (length(admin.cens)!=n) stop("censoring weights must have length equal to nrow in data\n");  
        Gcx <- admin.cens
        Gctimes <- rep(1,length(times)); 
    }

   if (left==1 & is.null(trunc.p) & is.null(cens.weights))  {  ## {{{ 
	 ### geskus weights: from mstate crprep 
	 stop("For left-truncated data call prep.comp.risk\n call with weights and cens.weights\n"); 
         n=length(time2)
         prec.factor <- 100
         prec <- .Machine$double.eps * prec.factor
         surv.trunc <- survfit(Surv(-time2,-entrytime+prec,rep(1,n)) ~ 1)
         trunc.dist <- summary(surv.trunc)
         trunc.dist$time <- rev(-trunc.dist$time)
         trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
         Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),time2)
         Lw <- Lfit[,2]
###	 weights <- 1/Lw
	 weights <- 1/((Lw)*Gcx); 
	 weights[delta==cens.code] <- 0
	 Gcx <- rep(1,n)
   } ## }}} 
   if (is.null(trunc.p)) trunc.p <- rep(1,n);  
   if (length(trunc.p)!=n) stop("truncation weights must have same length as data\n"); 
## }}}


## {{{ setting up more variables

  if (resample.iid == 1) {
    biid <- double(ntimes* antclust * px);
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  ps<-px; 
  betaS<-rep(0,ps); 

  ## possible starting value for nonparametric components
  if (is.null(est)) { 
	  est<-matrix(0.0+0.1,ntimes,px+1); 
          est[,1] <- times; 
  }  else {
      est <- as.matrix(est); 
  }
  if (nrow(est)!=length(times)) est <- Cpred(est,times); 

  hess<-matrix(0,ps,ps); var<-score<-matrix(0,ntimes,ps+1); 
  if (sum(gamma)==0) gamma<-rep(0,pg); gamma2<-rep(0,ps); 
  test<-matrix(0,antsim,3*ps); testOBS<-rep(0,3*ps); unifCI<-c();
  testval<-c(); rani<--round(runif(1)*10000); 
  Ut<-matrix(0,ntimes,ps+1); simUt<-matrix(0,ntimes,50*ps);
  var.gamma<-matrix(0,pg,pg); 
  pred.covs.sem<-0

  if (is.null(time.pow)==TRUE & model=="prop" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="fg" )       time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="additive")  time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif2" )    time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic" ) time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic2" )time.pow<-rep(0,pg); 
  if (length(time.pow)!=pg) time.pow <- rep(time.pow[1],pg); 

  if (is.null(time.pow.test)==TRUE & model=="prop" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="fg" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="additive")  time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif" )    time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif2" )   time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic" ) time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic2" ) time.pow.test<-rep(0,px); 
  if (length(time.pow.test)!=px) time.pow.test <- rep(time.pow.test[1],px); 

  if (ntimes>1) silent <- c(silent,rep(0,ntimes-1))
  ## }}}

###    print(Gctimes); 

  ###  dyn.load("comprisk.so"0
  ssf <- step;  ## takes step size over 
  out<-.C("itfit", ## {{{
          as.double(times),as.integer(ntimes),as.double(eventtime),
          as.integer(cens.code), as.integer(status),as.double(Gcx),
          as.double(X),as.integer(n),as.integer(px),
          as.integer(Nit), as.double(betaS), as.double(score),
          as.double(hess), as.double(est), as.double(var),
          as.integer(sim),as.integer(antsim),as.integer(rani),
          as.double(test), as.double(testOBS), as.double(Ut),
          as.double(simUt),as.integer(weighted),as.double(gamma),
          as.double(var.gamma),as.integer(fixed),as.double(Z),
          as.integer(pg),as.integer(trans),as.double(gamma2),
          as.integer(cause),as.integer(line),as.integer(detail),
          as.double(biid),as.double(gamiid),as.integer(resample.iid),
          as.double(time.pow),as.integer(clusters),as.integer(antclust),
          as.double(time.pow.test),as.integer(silent), as.double(conv),
	  as.double(weights),as.double(entrytime),as.double(trunc.p),
	  as.integer(estimator),as.integer(fix.gamma), as.integer(stratum),
	  as.integer(ordertime-1),as.integer(conservative), as.double(ssf), 
	  as.double(Gctimes),as.double(rep(0,pg)),as.double(matrix(0,pg,pg)),
	  as.integer(monotone),PACKAGE="timereg") ## }}}
 
 ## {{{ handling output
  ssf <- out[[51]]; 
  gamma<-matrix(out[[24]],pg,1); 
  var.gamma<-matrix(out[[25]],pg,pg); 
  Dscore.gamma<-matrix(out[[54]],pg,pg); 
  gamma2<-matrix(out[[30]],ps,1); 
  rownames(gamma2)<-covnamesX; 
  
  conv <- list(convp=out[[41]],convd=out[[42]]);
    
  if (fixed==0) gamma<-NULL; 

  if (resample.iid==1)  {
    biid<-matrix(out[[34]],ntimes,antclust*px);
    if (fixed==1) gamiid<-matrix(out[[35]],antclust,pg) else gamiid<-NULL; 
    B.iid<-list();
    for (i in (0:(antclust-1))*px) {
    B.iid[[i/px+1]]<-matrix(biid[,i+(1:px)],ncol=px);
      colnames(B.iid[[i/px+1]])<-covnamesX; }
    if (fixed==1) colnames(gamiid)<-covnamesZ
  } else B.iid<-gamiid<-NULL;

  if (sim==1) {
    simUt<-matrix(out[[22]],ntimes,50*ps); UIt<-list();
    for (i in (0:49)*ps) UIt[[i/ps+1]]<-as.matrix(simUt[,i+(1:ps)]);
    Ut<-matrix(out[[21]],ntimes,ps+1);
    test<-matrix(out[[19]],antsim,3*ps); testOBS<-out[[20]];
    supUtOBS<-apply(abs(as.matrix(Ut[,-1])),2,max);
    p<-ps
    for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
    for (i in 1:p) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
    pval.testBeq0<-as.vector(testval[1:p]);
    pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
    pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
    obs.testBeq0<-as.vector(testOBS[1:p]);
    obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
    obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
    sim.testBeq0<-as.matrix(test[,1:p]);
    sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
    sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
  } else {test<-unifCI<-Ut<-UIt<-pval.testBeq0<-pval.testBeqC<-obs.testBeq0<-
          obs.testBeqC<- sim.testBeq0<-sim.testBeqC<-
          sim.testBeqC.is<- pval.testBeqC.is<-
          obs.testBeqC.is<-NULL;
  }

    est<-matrix(out[[14]],ntimes,ps+1);
    ## in case of no convergence set estimates to NA 
    est[conv$convp>0,-1] <- NA
    score<-matrix(out[[12]],ntimes,ps+1); 
    gamscore <- matrix(out[[53]],pg,1)
    scores <- list(score=score,gamscore=gamscore)
    var<-matrix(out[[15]],ntimes,ps+1);
    ## in case of no convergence set var to NA 
    var[conv$convp>0,-1] <- NA
    colnames(var)<-colnames(est)<-c("time",covnamesX); 

  if (sim>=1) {
    colnames(Ut)<- c("time",covnamesX)
    names(unifCI)<-names(pval.testBeq0)<- names(pval.testBeqC)<- 
    names(pval.testBeqC.is)<- names(obs.testBeq0)<- names(obs.testBeqC)<- 
    names(obs.testBeqC.is)<- colnames(sim.testBeq0)<- colnames(sim.testBeqC)<- 
    colnames(sim.testBeqC.is)<- covnamesX;
  }

    if (fixed==1) { rownames(gamma)<-c(covnamesZ);
                    colnames(var.gamma)<- rownames(var.gamma)<-c(covnamesZ); }

    colnames(score)<-c("time",covnamesX);
    if (is.na(sum(score))==TRUE) score<-NA  else 
    if (sum(score[,-1])<0.00001) score<-sum(score[,-1]); 
  ud<-list(cum=est,var.cum=var,gamma=gamma,score=score,
           gamma2=gamma2,var.gamma=var.gamma,robvar.gamma=var.gamma,
           pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
           obs.testBeq0=obs.testBeq0,
           obs.testBeqC.is=obs.testBeqC.is,
           obs.testBeqC=obs.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
           conf.band=unifCI,B.iid=B.iid,gamma.iid=gamiid,ss=ssf,
           test.procBeqC=Ut,sim.test.procBeqC=UIt,conv=conv,
	   weights=weights,cens.weights=cens.weights,scores=scores,Dscore.gamma=Dscore.gamma,step=step)

    ud$call<-call; 
    ud$model<-model; 
    ud$n<-n; 
    ud$clusters <- clusters
    ud$formula<-formula;
    ud$response <- event.history
    ud$cause <- status
    class(ud)<-"comprisk"; 
    attr(ud, "Call") <- call
    attr(ud, "Formula") <- formula
    attr(ud, "time.pow") <- time.pow
    attr(ud, "causeS") <- causeS
    attr(ud, "cause") <- status
    attr(ud, "cluster.call") <- cluster.call
    attr(ud, "coarse.clust") <- coarse.clust
    attr(ud, "max.clust") <- max.clust
    attr(ud, "clusters") <- clusters
    attr(ud, "cens.code") <- cens.code
    attr(ud, "times") <- times
    return(ud);  ## }}}
} ## }}}

print.comprisk <- function (x,...) { ## {{{
  object <- x; rm(x);
  if (!inherits(object, 'comprisk')) stop ("Must be an comprisk object")
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:
  causeS <- attr(object,"causeS")
  print(causeS)
  cat(paste("\nAnalysed cause:",causeS,"\n"))      
  cat(paste("\nLink _function:",object$model,"\n\n"))
  cat(" Nonparametric terms : ");
  cat(colnames(object$cum)[-1]); cat("   \n");  
  if (semi) {
      cat(" Parametric terms :  ");
      cat(rownames(object$gamma)); 
      cat("   \n");
  } 

  if (object$conv$convd>=1) {
      if (all(object$conv$convp==1)){
          if (NROW(object$cum)>1){
              cat("\nWarning: problem with convergence at all time points\n")
          } else{
              cat("\nWarning: problem with convergence at the evaluation time.\n")
          }
      }else{
          cat("Warning: problem with convergence at the following time points:\n")
          cat(object$cum[object$conv$convp>0,1])
          cat("\nYou may try to readjust analyses by removing these time points\n")
      }
  }
  cat("   \n");  
} ## }}}

coef.comprisk <- function(object, digits=3,...) { ## {{{
   coefBase(object,digits=digits)
} ## }}}

summary.comprisk <- function (object,digits = 3,...) {  ## {{{
  if (!inherits(object, 'comprisk')) stop ("Must be a comprisk object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Competing risks Model \n\n")
  
  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 

  if (semi) { 
         if (sum(abs(object$score)>0.000001)) 
         cat("Did not converge, allow more iterations\n\n"); 

	 cat("Parametric terms : \n"); 
	 prmatrix(coef(object,digits=digits))
	 cat("   \n"); 
  }

  if (object$conv$convd>=1) {
       cat("WARNING problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n\n") }

} ## }}}

plot.comprisk <-  function (x, pointwise.ci=1, hw.ci=0,
                            sim.ci=0, specific.comps=FALSE,level=0.05, start.time = 0,
                            stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
                            ylab ="Coefficients",score=FALSE,...){
## {{{
  object <- x; rm(x);

  if (!inherits(object,'comprisk') ){
    stop ("Must be output from comp.risk function")
  }

  if (score==FALSE) {
    B<-object$cum;
    V<-object$var.cum;
    p<-dim(B)[[2]]; 

    if (sum(specific.comps)==FALSE){
      comp<-2:p
    } else {
      comp<-specific.comps+1
    }
    if (stop.time==0) {
      stop.time<-max(B[,1]);
    }

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,];
    V<-V[med,]; 

    c.alpha<- qnorm(1-level/2)
    for (v in comp) { 
      c.alpha<- qnorm(1-level/2)
      est<-B[,v];
      ul<-B[,v]+c.alpha*V[,v]^.5;
      nl<-B[,v]-c.alpha*V[,v]^.5;
      if (add.to.plot==FALSE) {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab,...) 
        if (mains==TRUE) title(main=colnames(B)[v]);
      } else {
        lines(B[,1],est,type="s");
      }
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s");
        lines(B[,1],nl,lty=pointwise.ci,type="s");
      }
      if (hw.ci>=1) {
        if (level!=0.05){
          cat("Hall-Wellner bands only 95 % \n");
        }
        tau<-length(B[,1])
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s");
      }
      if (sim.ci>=1) {
        if (is.null(object$conf.band)==TRUE){
          cat("Uniform simulation based bands only computed for n.sim> 0\n")
        }
        if (level!=0.05){
          c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        } else {
          c.alpha<-object$conf.band[v-1];
        }
        nl<-B[,v]-c.alpha*V[,v]^.5;
        ul<-B[,v]+c.alpha*V[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s");
      }
      abline(h = 0)
    }
  } else {
    # plot score proces
    if (is.null(object$pval.testBeqC)==TRUE) {
      cat("Simulations not done \n"); 
      cat("To construct p-values and score processes under null n.sim>0 \n"); 
    } else {
      if (ylab=="Cumulative regression function"){ 
        ylab<-"Test process";
      }
      dim1<-ncol(object$test.procBeqC)
      if (sum(specific.comps)==FALSE){
        comp<-2:dim1
      } else {
        comp<-specific.comps+1
      }

      for (i in comp){
          ranyl<-range(object$test.procBeqC[,i]);
          for (j in 1:50){
            ranyl<-range(c(ranyl,(object$sim.test.procBeqC[[j]])[,i-1]));
          }
          mr<-max(abs(ranyl));

          plot(object$test.procBeqC[,1],
               object$test.procBeqC[,i],
               ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab=ylab,type="s",...)
          if (mains==TRUE){
            title(main=colnames(object$test.procBeqC)[i]);
          }
          for (j in 1:50){
            lines(object$test.procBeqC[,1],
                  as.matrix(object$sim.test.procBeqC[[j]])[,i-1],col="grey",lwd=1,lty=1,type="s")
          }
          lines(object$test.procBeqC[,1],object$test.procBeqC[,i],lwd=2,type="s")
        }
    }
  }
} ## }}}

prep.comp.risk <- function(data,times=NULL,entrytime=NULL,
			   time="time",cause="cause",cname="cweight",tname="tweight",
			   strata=NULL,nocens.out=TRUE,cens.formula=NULL,cens.code=0,
			   prec.factor=100)
{ ## {{{ 
## {{{  geskus weights, up to min(T_i,max(times))
   if (is.null(times)) times <- max(data[,time])
   if (is.null(entrytime)) entrytime <- rep(0,nrow(data)) else entrytime <- data[,entrytime]
   mtt <- max(times)
   prec.factor <- 100
   prec <- .Machine$double.eps * prec.factor
   trunc.model <- cens.model <- NULL ## output of Cox models for entry cens

   if (is.null(cens.formula)) { 
   if (is.null(strata)) { ## {{{ 
	   if (!is.null(entrytime)) {
	   surv.trunc <- 
	   survfit(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ 1) 
	   trunc.dist <- summary(surv.trunc)
	   trunc.dist$time <- rev(-trunc.dist$time)
	   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
	   Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),data[,time])
	   Lw <- Lfit[,2]
	   } else Lw <- 1
	   ud.cens<- survfit(Surv(entrytime,data[,time],data[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,data[,time]),strict=TRUE)[,2];
           weights <- 1/(Lw*Gcx); 
	   cweights <-  Lw; 
	   tweights <-  Gcx; 
   ### ## }}} 
   } else { ## {{{ 
	   ### compute for each strata and combine 
	  vstrata <- as.numeric(data[,strata])
          weights <- rep(1,nrow(data))
          cweights <- rep(1,nrow(data))
          tweights <- rep(1,nrow(data))
	  for (i in unique(vstrata)) { ## {{{ for each strata
	       who <- (vstrata == i)
	       if (sum(who) <= 1) stop(paste("strata",i,"less than 1 observation\n")); 
	   datas <- subset(data,who)
	   if (!is.null(entrytime)) {
		   entrytimes <- entrytime[who]
		   surv.trunc <- 
		   survfit(Surv(-datas[,time],-entrytimes+prec,rep(1,nrow(datas))) ~ +1) 
		   trunc.dist <- summary(surv.trunc)
		   trunc.dist$time <- rev(-trunc.dist$time)
		   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
		   Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),datas[,time])
		   Lw <- Lfit[,2]
	   } else Lw <- 1
	   ud.cens<- survfit(Surv(entrytimes,datas[,time],datas[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,datas[,time]),strict=TRUE)[,2];
	   weights[who]<-  1/(Lw*Gcx); 
	   cweights[who]<-  Lw; 
	   tweights[who]<-  Gcx; 
          } ## }}} 
   } ## }}} 
   } else { ### cens.formula Cox models  ## {{{
        X <- model.matrix(cens.formula,data=data)[,-1,drop=FALSE]; 

	if (!is.null(entrytime)) {
		trunc.model <- coxph(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ X) 
		baseout <- basehaz(trunc.model,centered=FALSE); 
		baseout <- cbind(rev(-baseout$time),rev(baseout$hazard))
	###
		Lfit <-Cpred(baseout,data[,time])[,-1]
		RR<-exp(as.matrix(X) %*% coef(trunc.model))
		Lfit<-exp(-Lfit*RR)
		Lw <- Lfit
	   } else Lw <- 1
###
	cens.model <- coxph(Surv(entrytime,data[,time],data[,cause]==0)~+X) 
        baseout <- basehaz(cens.model,centered=FALSE); 
	baseout <- cbind(baseout$time,baseout$hazard)
	Gfit<-Cpred(baseout,pmin(mtt,data[,time]),strict=TRUE)[,2];
	RR<-exp(as.matrix(X) %*% coef(cens.model))
	Gfit<-exp(-Gfit*RR)
        weights <- 1/(Lw*Gfit); 
        cweights <- Gfit
        tweights <- Lw
   } ## }}} 
   data[,cname] <- cweights
   data[,tname] <- tweights

   if (!is.null(entrytime)) {
   mint <- min(tweights); maxt <- min(tweights) 
   if (mint<0 | mint>1) warning("min(truncation weights) strange, maybe prec.factor should be different\n")
   if (maxt<0 | maxt>1) warning("max(truncation weights) strange, maybe prec.factor should be different\n")
   }

   if ("weights" %in% names(data)) {
       warning("Weights in variable 'weights_' \n")
       wname<- "weights_"
       data[,wname] <- weights
   } else data[,"weights"] <- weights
###
   if ("cw" %in% names(data)) {
     warning("cw weights in variable 'cw_' \n")
     cwname<- "cw_"
     data[,cwname] <- 1
   } else data[,"cw"] <- 1
###
   if (nocens.out) {
     med <- ((data[,time]>mtt & data[,cause]==cens.code)) | (data[,cause]!=cens.code)
     data <- data[med,]
   } 

   attr(data,"trunc.model") <- trunc.model
   attr(data,"cens.model") <- cens.model 
## }}} 
   return(data)
} ## }}} 

pred.stratKM <- function(data,entrytime=NULL,time="time",cause="cause",strata="strata",event.code=0)
{ ## {{{ 

     if (is.numeric(time)) time <- time else {
	     if (!is.null(data)) time <- data[,time] else stop("time not given\n"); 
     }
     if (is.numeric(cause)) cause <- cause else {
	     if (!is.null(data)) cause <- data[,cause] else stop("cause not given\n"); 
     }
     if (is.numeric(strata)) strata <- strata else {
	     if (!is.null(data)) strata <- data[,strata] else stop("strata not given\n"); 
     }
  if (is.null(entrytime)) entrytime <- rep(0,nrow(data)) else {
     if (is.numeric(entrytime)) entrytime <- entrytime else {
	     if (!is.null(data)) entrytime <- data[,entrytime] else stop("entrytime not given\n"); 
     }
  }
  vstrata <- as.numeric(strata)
  weights <- rep(1,length((data)))
  for (i in unique(vstrata)) { ## {{{ for each strata
	   who <- (vstrata == i)
	   if (sum(who) <= 1) stop(paste("strata",i,"less than 1 observation\n")); 
	   times <- time[who]
	   causes <- cause[who]
	   entrytimes <- entrytime[who]
	   ud.cens<- survfit(Surv(entrytimes,times,causes==event.code)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,times,strict=TRUE)[,2];
	   weights[who]<-  Gcx; 
   } ## }}} 
   return(weights); 
} ## }}} 

