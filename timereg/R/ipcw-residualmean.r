res.mean<-function(formula,data=sys.parent(),cause=1,restricted=NULL,times=NULL,Nit=50,
clusters=NULL,gamma=0,n.sim=0,weighted=0,model="additive",detail=0,interval=0.01,resample.iid=1,
cens.model="KM",cens.formula=NULL,time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6,estimator=1,cens.weights=NULL,
conservative=1,weights=NULL){
## {{{
# restricted residual mean life models, using IPCW 
# two models, additive and proportional 
# trans=1 E(min(Y,tau) - t | Y>= t) = ( x' b(b)+ z' gam (tau-t)) 
# trans=2 E(min(Y,tau) - t | Y>= t) = exp( x' b(b)+ z' gam (tau-t)) 
# unrestricted residual mean life models, using IPCW 
# trans=1 E(Y - t | Y>= t) = ( x' b(b)+ z' gam ) 
# trans=2 E(Y - t | Y>= t) = exp( x' b(b)+ z' gam ) 
  if (model=="additive") trans<-1; 
  if (model=="prop")     trans<-2; 
  if (model=="additive") Nit=2;  ### one for estimation and one for residuals 
# estimator =1   E(min(Y,tau)-t | Y>=t)  or  E(Y - t | Y>=t) 
# estimator =2   Year lost due to cause 1 up to time tau given event
# estimator =2   E( tau - min(T,tau) | T <= tau, epsilon=j)
# estimator =3   PKA Years lost due to cause 1 up to time tau
# estimator =3   E( tau - min(T_j,tau) | T>t ) = \int_t^tau (1-F_j(s)) ds / S(t)
# estimator =3   E( tau - min(T,tau) | T <= tau, epsilon=j, T>t) * F_j(tau)
  cens.tau <- 1

  line <- 0
  m<-match.call(expand.dots = FALSE);
  m$gamma<-m$times<-m$cause<-m$Nit<-m$weighted<-m$n.sim<-
###	   m$cens.tau <- 
           m$model<- m$detail<- m$cens.model<-m$time.pow<-m$silent<- 
           m$interval<- m$clusters<-m$resample.iid<-m$restricted <- m$weights <- 
           m$time.pow.test<-m$conv<-m$estimator <- m$cens.weights <- 
	   m$conservative <- m$cens.formula <- NULL
  special <- c("const","cluster")
  if (missing(data)) {
    Terms <- terms(formula, special)
  }  else {
    Terms <- terms(formula, special, data = data)
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")

   ## {{{ Event stuff
    cens.code <- attr(Y,"cens.code")
    if (ncol(Y)==3) stop("Left-truncation, through weights \n"); 
    time2 <- eventtime <- Y[,1]
    status <- delta  <- Y[,2]
    event <- (status==cause)
    entrytime <- rep(0,length(time2))
    if (sum(event)==0) stop("No events of interest in data\n"); 
    ## }}} 

  if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;

  if(is.null(clusters)){ clusters <- des$clusters}

  if(is.null(clusters)){
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    clusters <- as.integer(factor(clusters))-1
    antclust <- length(unique(clusters))
  }
    
  pxz <-px+pz;
  ### always survival case 
  status; ### cause==1 and cause==0 censoring
  if (is.null(restricted)) tau <- max(time2)+0.1 else tau <- restricted

  if (is.null(times)) {
      times<-sort(unique(time2[status==cause])); 
      if (tau>0) times <- times[times<tau]; 
      ###times<-times[-c(1:5)];
  } else times <- sort(times); 

  n<-nrow(X); ntimes<-length(times);
  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
  ### non-cens  at restricted time
  delta<-(status!=cens.code)
  if (!is.null(restricted)) { 
     if (cens.tau==1) {
       time2tau <- pmin(time2,restricted)
       deltatau <- rep(1,n)
       deltatau[status==cens.code & time2<restricted] <- 0
      } else { time2tau <- time2; deltatau <- delta*1;}
  } else  { time2tau <- time2; deltatau <- delta*1;}

  if (is.null(weights)==TRUE) weights <- rep(1,n); 

###  cens model  ## {{{ 
  if (cens.model=="KM") { ## {{{ 
    ud.cens<-survfit(Surv(time2,status==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time2tau)[,2];
    Gtimes<-Cpred(Gfit,times)[,2];
    Gctimes<-Cpred(Gfit,times)[,2];
    Gctimes[Gctimes<=0] <- 1 ## }}} 
  } else if (cens.model=="cox") {  ## {{{ 
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-cox.aalen(Surv(time2,status==cens.code)~prop(XZ),n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time2tau)[,2];
    RR<-exp(XZ %*% ud.cens$gamma)
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2];
    Gctimes[Gctimes<=0] <- 1  ## }}} 
    } else if (cens.model=="aalen") {  ## {{{ 
    if (npar==TRUE) XZ <- X[,-1] else XZ <- cbind(X,Z)[,-1];
    ud.cens <- aalen(Surv(time2,status==cens.code)~XZ,n.sim=0,robust=0);
    Gcx <- Cpred(ud.cens$cum,time2tau)[,-1];
    Gcx1 <- Gcx
    XZ <- cbind(1,XZ); 
    Gcx <- exp(-apply(Gcx*XZ,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    Gctimes<-Cpred(Gfit,times)[,2];
    Gctimes[Gctimes<=0] <- 1 ## }}} 
    } else if (cens.model=="weights") { 
      if (length(weights)!=n) stop(paste("Weights length=",length(weights),"do not have length \n",length(time2),"problems with missing values?\n")); 
      Gcx <- cens.weights
      ord2 <- order(time2)
      Gctimes <- Cpred(cbind(time2[ord2],Gcx[ord2]),times)
      Gctimes[Gctimes<=0] <- 1
    } else { stop('Unknown censoring model') }
  ## }}} 

  if (resample.iid == 1) {
    biid <- double(ntimes* antclust * px);
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  ps<-px; 
  hess<-matrix(0,ps,ps); var<-score<-matrix(0,ntimes,ps+1); 
  if (sum(gamma)==0) gamma<-rep(0,pg); gamma2<-rep(0,ps); 
  test<-matrix(0,antsim,3*ps); testOBS<-rep(0,3*ps); unifCI<-c();
  testval<-c(); rani<--round(runif(1)*10000); 
  Ut<-matrix(0,ntimes,ps+1); simUt<-matrix(0,ntimes,50*ps);
  var.gamma<-matrix(0,pg,pg); 
  pred.covs.sem<-0

  if (is.null(time.pow)==TRUE & !is.null(restricted))       time.pow<-rep(1,pg); 
  if (is.null(time.pow.test)==TRUE & !is.null(restricted))  time.pow<-rep(1,pg); 
###  if (is.null(time.pow)==TRUE & model=="additive")  time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & is.null(restricted))       time.pow<-rep(0,pg); 
  if (is.null(time.pow.test)==TRUE & is.null(restricted))  time.pow<-rep(0,pg); 

  silent <- c(silent,rep(0,ntimes-1));
  if (!is.null(restricted)) time2  <- pmin(time2,restricted) 
###  print(restricted); print(time2); print(table(deltatau)); print(table(status)); 
###  print(table(status)); print(times); print(summary(Gctimes)); print(summary(Gcx)); 

  ### important to start in 0 for linear model 
  est<-matrix(0,ntimes,ps+1) 
  if (model!="additive") est[,2] <- log(mean(time2))
  betaS<-rep(0,ps);  
  if (model=="prop") betaS[1] <- log(mean(time2))

  ordertime <- order(eventtime); 

  out<-.C("resmean",
      as.double(times),as.integer(ntimes),as.double(time2),
      as.integer(deltatau), as.integer(status),as.double(Gcx),
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
      as.double(time.pow.test),as.integer(silent),
      as.double(conv),as.double(tau),as.integer(estimator),
      as.integer(cause),as.double(weights),as.double(Gctimes),
      as.integer(ordertime-1),as.integer(conservative),as.integer(cens.code),
      PACKAGE="timereg")

  gamma<-matrix(out[[24]],pg,1); 
  var.gamma<-matrix(out[[25]],pg,pg); 
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
  score<-matrix(out[[12]],ntimes,ps+1); 
  var<-matrix(out[[15]],ntimes,ps+1); 
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
           conf.band=unifCI,B.iid=B.iid,gamma.iid=gamiid,
           test.procBeqC=Ut,sim.test.procBeqC=UIt,conv=conv,
	   cens.weights=Gcx,time=time2,delta.tau=deltatau,time2tau=time2tau)

  ud$call<-call; ud$model<-model; ud$n<-n; 
  ud$formula<-formula; class(ud)<-"resmean"; 
  attr(ud, "Call") <- sys.call()
  attr(ud, "Formula") <- formula
  attr(ud, "time.pow") <- time.pow
  attr(ud, "cause") <- cause
  attr(ud, "restricted") <- restricted
  attr(ud, "times") <- times
  attr(ud, "model") <- model
  attr(ud, "estimator") <- estimator
  return(ud); 
} ## }}}

print.resmean <- function (x,...) { ## {{{
  object <- x; rm(x);
  if (!inherits(object, 'resmean')) stop ("Must be an resmean object")
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
   # We print information about object:  
  cat(paste("Residual mean model with",object$model,"\n"))
  if (!is.null(attr(object,"restricted")))
  cat(paste("Restricted at ",attr(object,"restricted","\n")))
  cat("\n")
  cat(" Nonparametric terms : ");
  cat(colnames(object$cum)[-1]); cat("   \n");  
  if (semi) {
    cat(" Parametric terms :  ");
    cat(rownames(object$gamma)); 
    cat("   \n");
  } 

  if (object$conv$convd>=1) {
       cat("Warning problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n") }
  cat("   \n");  
} ## }}}

coef.resmean <- function(object, digits=3,...) { ## {{{
   coefBase(object,digits=digits)
} ## }}}

summary.resmean <- function (object,digits = 3,ci=0, alpha=0.05,silent=0, ...) { ## {{{
  if (!inherits(object, 'resmean')) stop ("Must be a resmean object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  if (silent==0) {
  # We print information about object:  
  cat(paste("Residual mean model with",object$model,"\n"))
###  cat("Residual mean model \n\n")
  if (attr(object,"estimator")==1) cat("Residual mean life:")
  if (attr(object,"estimator")==2) cat("Years Lost to cause given event:")
  if (attr(object,"estimator")==3) cat("Years Lost to cause:")
  if (!is.null(attr(object,"restricted")))
  cat(paste("Restricted at ",attr(object,"restricted","\n")))
  cat("\n"); cat("\n")

  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 
  }

  if (semi) { if (silent==0) cat("Parametric terms : \n"); 
              out=coef.resmean(object); 
              out=signif(out,digits=digits)
	      if (silent==0) { print(out); cat("   \n"); }
  } else { if (silent==0) cat("Non-Parametric terms : \n"); 
          if (nrow(object$cum)==1) {
	     out=cbind(c(object$cum[,-1]),c(object$var.cum[,-1]^.5))
	     colnames(out) <- c("resmean","se")
             out=signif(out,digits=digits)
	  } else {
             se.cum <- object$var.cum[,-1]^.5
	     colnames(se.cum) <- paste("se",colnames(se.cum),paste="")
	     out=cbind(object$cum,se.cum)
             out=signif(head(out),digits=digits)
	  }
   if (silent==0)       cat("   \n"); 
  }
  if (ci==1) { out <- round(cbind(out,out[,1]+qnorm(alpha/2)*out[,2],out[,1]+qnorm(1-alpha/2)*out[,2]),digits=digits); 
               nn <- ncol(out); 
               colnames(out)[(nn-1):nn] <- c("lower","upper")
  }

  if (silent==0) {
  if (object$conv$convd>=1) {
       cat("WARNING problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n\n") 
  }

  cat("  Call: \n")
  dput(attr(object, "Call"))
  cat("\n")
  }

  out
} ## }}}

plot.resmean <-  function (x, pointwise.ci=1, hw.ci=0,
                            sim.ci=0, specific.comps=FALSE,level=0.05, start.time = 0,
                            stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
                            ylab ="Coefficients",score=FALSE,...){ ## {{{
  object <- x; rm(x);

  if (!inherits(object,'resmean') ){ stop ("Must be output from res.mean function") }

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
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="l",xlab=xlab,ylab=ylab) 
        if (mains==TRUE) title(main=colnames(B)[v]);
      } else {
        lines(B[,1],est,type="l");
      }
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="l");
        lines(B[,1],nl,lty=pointwise.ci,type="l");
      }
      if (hw.ci>=1) {
        if (level!=0.05){
          cat("Hall-Wellner bands only 95 % \n");
        }
        tau<-length(B[,1])
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="l"); 
        lines(B[,1],nl,lty=hw.ci,type="l");
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
        lines(B[,1],ul,lty=sim.ci,type="l"); 
        lines(B[,1],nl,lty=sim.ci,type="l");
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
               ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab=ylab,type="l")
          if (mains==TRUE){
            title(main=colnames(object$test.procBeqC)[i]);
          }
          for (j in 1:50){
            lines(object$test.procBeqC[,1],
                  as.matrix(object$sim.test.procBeqC[[j]])[,i-1],col="grey",lwd=1,lty=1,type="l")
          }
          lines(object$test.procBeqC[,1],object$test.procBeqC[,i],lwd=2,type="l")
        }
    }
  }
} ## }}}

