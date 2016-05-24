cum.residuals<-function(object,data=sys.parent(),modelmatrix=0,cum.resid=1,n.sim=500,
	weighted.test=0,max.point.func=50,weights=NULL)
{ ## {{{
## {{{ setting up
  start.design<-1; silent <- 1; offsets <- NULL; 
  if (!(class(object)!="aalen" | class(object)!="timecox" | class(object)!="cox.aalen" ))
    stop ("Must be output from aalen() timecox() or cox.aalen() functions\n") 
  if (class(object)=="timecox") if (object$method!="basic") 
    stop("Residuals available only for method=basic\n")
  if (class(object)=="timecox") if (is.null(object$gamma)==FALSE) 
    stop("Residuals available only for timecox model with no const terms\n")
  if (class(object)=="aalen") if (is.null(object$gamma)==FALSE) 
    stop("Residuals available only for Aalen model with no const terms\n")
  if (is.null(object$residuals$dM)==TRUE) stop("Residuals not computed, add option residuals=1\n");
  if (sum(modelmatrix)==0 && cum.resid==0) 
	  stop("No modelmatrix or continous covariates given to cumulate residuals\n"); 
  stratum <- attr(object,"stratum"); 
  rate.sim <- 1; 
  weights1 <- NULL
  if (class(object)=="cox.aalen") {
    dcum<-apply(as.matrix(object$cum[,-1]),2,diff); 
    beta<-object$gamma; 
    coxaalen<-1; 
    weights1 <- attr(object,"weights") 
    offsets <- attr(object,"offsets"); 
    rate.sim <- attr(object,"rate.sim"); 
  } else { dcum<-0; beta<-0; coxaalen<-0; pg<-0;Z<-0; }

  id<-attr(object,"id"); 
  cluster<-attr(object,"cluster"); 
  formula<-attr(object,"Formula"); 
  start.time<-attr(object,"start.time"); 
  pers<-unique(id); 
  antpers<-length(pers); 
  clust<-unique(cluster); 
  antclust<-length(clust); 

  if (class(object)=="cox.aalen") 
  ldata<-aalen.des(formula,data,model="cox.aalen") else ldata<-aalen.des(formula,data) 

  X<-ldata$X; covar<-X; px<-ldata$px; 

  time<-attr(object,"start"); 
  time2<-attr(object,"stop"); 
  if (sum(time)==0) type <- "right" else type <- "counting"
  status<-attr(object,"status");  
  if (is.null(weights)) weights <- rep(1,nrow(X));  
  if (is.null(weights1)) weights1 <- rep(1,nrow(X));  
  if (is.null(offsets)) offsets <- rep(0,nrow(X));  
  if (length(weights)!=nrow(X)) stop("Lengths of weights and data do not match\n"); 
  if (length(weights1)!=nrow(X)) stop("Lengths of weights from aalen/cox.aalen and data do not match\n"); 

  if (coxaalen==1) {
    Z<-ldata$Z;  covnamesZ<-ldata$covnamesZ; 
    pg<-ncol(Z); 
  } else Z <- 0
  Ntimes <- sum(status); 

  if (sum(modelmatrix)==0) {modelmatrix<-0;model<-0;pm<-1;}  else 
  {model<-1; modelmatrix<-as.matrix(modelmatrix); pm<-ncol(modelmatrix); 
   test<-matrix(0,n.sim,3*pm); testOBS<-rep(0,2*pm); 
   covnames<-colnames(modelmatrix); 
  }

###  print(id); print(cluster)

  times<-c(start.time,time2[status==1]); 
  times<-sort(times);
  antpers=length(unique(id)); 
  ntot<-nrow(X); 
  lmgresids<-length(object$residuals$time); 

###if ( type == "right" )  {  ## {{{
   ot<-order(-time2,status==1); # order in time, status=0 first for ties
   time2<-time2[ot]; 
   status<-status[ot]; 
   X<-as.matrix(X[ot,])
   if (coxaalen==1) Z<-as.matrix(Z[ot,])
   if (model==1) modelmatrix<-as.matrix(modelmatrix[ot,])
   start <- time[ot] ### fra call 
   stop<-time2;      ### fra call
   cluster<-cluster[ot]
   id<-id[ot];
   weightsmg <- weights[ot]
   weights <- weights1[ot]
   offsets <- offsets[ot]
   entry=rep(-1,ntot); 
###  } else {
###        eventtms <- c(time,time2)
###        status <- c(rep(0, ntot), status)
###        ix <- order(-eventtms,status==1)
###        etimes    <- eventtms[ix]  # Entry/exit times
###	status <- status[ix]
###        stop  <- etimes; 
###        start <- c(time,time)[ix]; 
###        tdiff    <- c(-diff(etimes),start.time) # Event time differences
###        entry  <- c(rep(c(1, -1), each = ntot))[ix]
###        X        <- as.matrix(X[rep(1:ntot, 2)[ix],])
###	if (coxaalen==1) Z <- as.matrix(Z[rep(1:ntot,2)[ix],])
###        if (model==1) modelmatrix<-as.matrix(modelmatrix[rep(1:ntot,2)[ix],])
###	id <- rep(id,2)[ix]
###	cluster <- rep(cluster,2)[ix]
###        weights <- rep(weights, 2)[ix]
###	offsets <- rep(offsets,2)[ix]
###    } ## }}}

  ntot <- nrow(X); 

  if (coxaalen==1) { gamma.iid<-object$gamma.iid 
                     covar<-cbind(X,Z);
                     ptot<-px+pg;
  } else { covar<-covar; ptot<-px; gamma.iid<-0; }

  covnames0<-colnames(covar); 
  if (is.null(covnames0)) covnames0  <- rep("",ptot); 
  antal<-0;  maxval<-0; intercept<-0; antal2<-0;  maxval2<-0; 
  xvals<-list();  ant<-rep(0,ptot);  
  for (i in 1:ptot) xvals[[i]]<-c(); 
  rani<-round(runif(1)*10000);
  keepcumz<-c(); k<-0
  for (j in 1:ptot) { 
    z<-unique(covar[,j]); z<-sort(z); 
    if (length(z)> max.point.func) z <- quantile(z,probs=seq(0,1,length=max.point.func))
    antal<-antal+1; ant[j]<-length(z); 
    if (ant[j]>2) { k<-k+1; keepcumz<-c(keepcumz,j); xvals[[k]]<-z;
        maxval<-max(maxval,length(z)); }
  }

  if (sum(keepcumz)==0 && cum.resid==1) 
    stop(" No continous covariates given to cumulate residuals \n"); 

  pcumz<-length(keepcumz); 
  uni.test<-matrix(0,n.sim,pcumz); 
  univar.proc<-matrix(0,maxval,pcumz);
  robvarcumz<-matrix(0,maxval,pcumz);
  sim.univar.proc<-matrix(0,maxval,50*pcumz); 
  time.proc<-matrix(0,lmgresids,2); 
  sim.time.proc<-matrix(0,lmgresids,50); 
  simcumz<-matrix(0,n.sim,pcumz); 
  xval<-matrix(0,maxval,pcumz); 
  k<-1; 
  for (i in keepcumz) {xval[1:ant[i],k]<-xvals[[k]]; k<-k+1;}

  # testOBS size and location of supremum, test simulated sup's
  unitime.test<-time.test<-mult.test<-multtime.test<- matrix(0,n.sim,2*pcumz); 
  unitime.testOBS<-uni.testOBS<-time.testOBS<-mult.testOBS<-
  multtime.testOBS<-rep(0,pcumz); 

  inXorZ<-rep(0,pcumz); 

  if (pcumz>0) {
  antX<-(sum(ant[1:px]>2)) 
    if (antX>0) {inX<-(1:px)[(ant[1:px]>2)]; inXorZ[1:antX]<-1} else {inX<-c();} 
    if (coxaalen==1) {inZ<-(1:pg)[(ant[(px+1):(px+pg)]>2)]; }  else inZ<-c()
    inXZ<-c(inX,inZ); 
    inXZ<-inXZ-1
  } else {inXZ<-0; inXorZ<-0} 
  ant<-ant[keepcumz]; 

  Ut<- cummgt<- robvarcum <- matrix(0,lmgresids,pm+1); 
  simUt<-matrix(0,lmgresids,pm*50); test<-matrix(0,n.sim,3*pm); 
  testOBS <- rep(0,3*pm); 
  ## }}}

###dyn.load("mgresid.so");
  dNit <- 0
  if (coxaalen==1) dNit <- object$residuals$dNit

  mgout<- .C("mgresid", ## {{{
     as.double(X),as.integer(ntot),as.integer(px), 
     as.integer(antpers),as.double(start),as.double(stop),
     as.integer(status),as.integer(id),as.double(object$residuals$time),
     as.integer(lmgresids),as.double(object$residuals$dM),as.integer(n.sim),
     as.double(xval), as.integer(ant), as.double(univar.proc),
     as.double(time.proc),as.double(sim.univar.proc), as.double(sim.time.proc),
     as.double(uni.test),as.double(uni.testOBS), as.double(time.test),
     as.double(time.testOBS),as.double(unitime.test), as.double(unitime.testOBS),
     as.double(modelmatrix),as.integer(model), as.integer(pm),
     as.double(cummgt),as.double(dNit), as.double(robvarcum),  # 10 
     as.double(testOBS),as.double(test), as.double(simUt),
     as.double(Ut),as.integer(cum.resid), as.integer(maxval),
     as.integer(start.design),as.integer(coxaalen), as.double(dcum),
     as.double(beta),as.double(Z), as.integer(pg),
     as.double(gamma.iid),as.integer(cluster), as.integer(antclust), 
     as.double(robvarcumz), as.double(simcumz), as.integer(inXZ), 
     as.integer(inXorZ),as.integer(pcumz), as.integer(entry),
     as.integer(stratum),as.integer(silent),as.double(weights1),
     as.double(offsets),as.integer(rate.sim),as.double(weights),
     as.integer(weighted.test)) ### , PACKAGE="timereg")
## }}}

## {{{ handling output from C
  if (model==1) {
    cum<-matrix(mgout[[28]],lmgresids,pm+1); 
    robvar.cum<-matrix(mgout[[30]],lmgresids,pm+1);
    var.cum<-robvar.cum; Ut<-matrix(mgout[[34]],lmgresids,pm+1);

    colnames(Ut)<-colnames(cum)<-colnames(var.cum)<- colnames(robvar.cum)<- c("time",covnames)
    test.procBeq0<-Ut; 

      simUt<-matrix(mgout[[33]],lmgresids,50*pm);
      UIt<-list(); 
      for (i in (0:49)*pm) {
	UIt[[i/pm+1]]<-as.matrix(simUt[,i+(1:pm)]);
      }
      testOBS<-mgout[[31]];
      test<-matrix(mgout[[32]],n.sim,3*pm);
      testval<-c(); unifCI<-c(); 
      for (i in 1:(2*pm)) testval<-c(testval,pval(test[,i],testOBS[i]))
      for (i in 1:pm) unifCI<-as.vector(c(unifCI,percen(test[,2*pm+i],0.95)));
      obs.testBeq0<-as.vector(testOBS[1:pm]);
      obs.testBeq0.is<-as.vector(testOBS[(pm+1):(2*pm)]);
      pval.testBeq0<-as.vector(testval[1:pm]);
      pval.testBeq0.is<-as.vector(testval[(pm+1):(2*pm)]);
      sim.testBeq0<-test[,(2*pm+1):(3*pm)]; 

      sim.test.procBeq0<-UIt; 
      names(unifCI)<- names(pval.testBeq0)<- names(obs.testBeq0)<- 
        names(pval.testBeq0.is)<- names(obs.testBeq0.is)<- covnames
  } else {
    cum<-robvar.cum<-test<-unifCI<-Ut<-UIt<-pval.testBeq0<-
      pval.testBeq0.is<-obs.testBeq0<-obs.testBeq0.is<-sim.testBeq0<-NULL; 
  } 

  if (cum.resid>=1) {
    univar.p<-matrix(mgout[[15]],maxval,pcumz)
    robvarcumz<-matrix(mgout[[46]],maxval,pcumz)
    simcumz<-matrix(mgout[[47]],n.sim,pcumz)
    univar.proc<-list(); 
    for (i in 1:pcumz) {
      univar.proc[[i]]<-cbind(xvals[[i]],univar.p[1:ant[i],i]); 
      colnames(univar.proc[[i]])<-c(covnames0[keepcumz[i]],"cum. martingale residual"); 
    }
    Uiz<-matrix(mgout[[17]],maxval,50*pcumz);
    UIz<-list(); 
    k<-1; 
    for (i in 1:pcumz) { UIz[[i]]<-matrix(Uiz[1:ant[i],i+(0:49)*pcumz],ncol=50); k<-k+1;}
    uni.test<-matrix(mgout[[19]],n.sim,pcumz)
    uni.test<-as.matrix(uni.test); 
    uni.testOBS<-mgout[[24]][1:pcumz]; 
    testval<-c(); 

    for (i in 1:pcumz)  testval<-c(testval,pval(uni.test[,i],uni.testOBS[i])) 
    unifCIz<-c()
    for (i in 1:pcumz)  unifCIz<-c(unifCIz,percen(simcumz[,i],0.95))

    uni.pval<-testval
    names(uni.testOBS)<-names(uni.pval)<-colnames(uni.test)<-covnames0[keepcumz]; 
    #uni.testOBS<-uni.testOBS[keepcumz]; uni.pval<-uni.pval[keepcumz]
    #uni.test<-uni.test[,keepcumz];
  } else { unifCIz<-uni.testOBS<-uni.pval<-proc.cumz<-UIz<-
           unitime.pval<-unitime.testOBS<-NULL;}

## }}}

  ud<-list(cum=cum,robvar.cum=robvar.cum,robvar.cumz=robvarcumz,
      pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0,
      pval.testBeq0.is=pval.testBeq0.is, obs.testBeq0.is=obs.testBeq0.is,
      sim.testBeq0=sim.testBeq0, procBeq0=Ut,sim.test.procBeq0=UIt,
      conf.band=unifCI, conf.band.cumz=unifCIz,
      obs.test=uni.testOBS,pval.test=uni.pval, sim.test=uni.test,
      proc.cumz=univar.proc,sim.test.proccumz=UIz)

  attr(ud,"Call")<-call; 
  class(ud)<-"cum.residuals"
  return(ud); 
} ## }}}

"print.cum.residuals"<- function (x,...)
{ ## {{{
    object <- x; rm(x);
	if (!inherits(object, 'cum.residuals')) stop ("Must be an MG resid object")

	cat("  Call: \n")
	dput(attr(object, "Call"))
	cat("\n")
} ## }}}

"plot.cum.residuals" <- function (x,pointwise.ci=1,hw.ci=0,sim.ci=0,
	robust=1, specific.comps=FALSE,level=0.05, start.time = 0, 
	stop.time = 0, add.to.plot=FALSE, mains=TRUE, main=NULL,
	xlab=NULL,ylab ="Cumulative MG-residuals",ylim=NULL,
	score=0,conf.band=FALSE,...) 
{## {{{
  object <- x; rm(x);
  if (!inherits(object,'cum.residuals') ) stop ("Must be output from cum.residuals()") 

  if (score <2) {
    B<-object$cum;
    if (sum(B)==0)  {
      stop("To compute cumulative residuals provide model matrix \n");
    } 
  }
  if (score==2)  ## {{{
  {
    if (sum(object$obs.test)==0)
      stop("To plot cumulative residuals vs. covariates, cum.resid=1"); 
  }

  if (score==0)  ## {{{ 
  {
    B<-object$cum; 
    if (robust>=1) V<-object$robvar.cum else V<-object$var.cum
    p <- ncol(B); 
    if (!is.null(main)) { if (length(main)!=p) main <- rep(main,length(comp)); mains <- FALSE; } 
    if (!is.null(xlab)) { if (length(xlab)!=p) xlab  <- rep(xlab,length(comp)); } 

    if (specific.comps==FALSE) comp<-(2:p) else comp<-specific.comps+1
    if (stop.time==0) stop.time<-max(B[,1]);

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
    V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
    Vrob<-object$robvar.cum; 
    Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

    c.alpha<- qnorm(1-level/2)
    for (v in comp)  ## {{{ 
    { 
      c.alpha<- qnorm(1-level/2)
      est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
      if (add.to.plot==FALSE) 
      {
       if (is.null(xlab)) xlabl <- "Time" else xlabl <- xlab[v]

       if (is.null(ylim)) plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlabl,ylab=ylab,...) 
       else plot(B[,1],est,ylim=ylim,type="s",xlab=xlabl,ylab=ylab) 

       if (!is.null(main)) title(main=main[i]); 
       if (mains==TRUE) title(main=colnames(B)[v]); 
      } else lines(B[,1],est,type="s"); 

      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s");
        lines(B[,1],nl,lty=pointwise.ci,type="s"); 
      }
      if (robust>=1) {
        lines(B[,1],ul,lty=robust,type="s"); 
        lines(B[,1],nl,lty=robust,type="s"); 
      }
      if (hw.ci>=1) {
        if (level!=0.05) cat("Hall-Wellner bands only 95 % \n");
        tau<-length(B[,1])
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s"); 
      }
      if (sim.ci>=1) {
        if (level!=0.05) c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        else c.alpha<-object$conf.band[v-1];
        nl<-B[,v]-c.alpha*Vrob[,v]^.5; ul<-B[,v]+c.alpha*Vrob[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s"); 
      }
      abline(h=0)
    } ## }}}

  }  ## }}} 
  else if (score==1)  ## {{{
  {
      dim1<-ncol(object$procBeq0)
      if (sum(specific.comps)==FALSE) comp<-2:dim1 else comp<-specific.comps+1

      if (!is.null(main)) { if (length(main)==1) main <- rep(main,length(comp)); mains <- FALSE; } 
      if (!is.null(xlab)) { if (length(xlab)==1) xlab <- c("time",rep(xlab,length(comp))); } 

      for (i in comp) ## {{{ 
        {
          ranyl<-range(object$procBeq0[,i]);
          for (j in 1:50) ranyl<-range(c(ranyl,(object$sim.test.procBeq0[[j]])[,i-1]));
          mr<-max(abs(ranyl));

          if (add.to.plot==FALSE)
	  {
             if (is.null(xlab)) xlabl <- "Time" else xlabl <- xlab[i]

	     if (!is.null(ylim)) 
             plot(object$procBeq0[,1],object$procBeq0[,i],type="s",
              ylim=ylim,lwd=2,xlab=xlabl,ylab=ylab,...)
	     else 
             plot(object$procBeq0[,1],object$procBeq0[,i],type="s",
                ylim=c(-mr,mr),lwd=2,xlab=xlabl,ylab=ylab,...)

            if (!is.null(main)) title(main=main[i]); 
            if (mains==TRUE) title(main=colnames(B)[i]); 
	  }
          else lines(object$procBeq0[,1],object$procBeq0[,i],type="s")


          for (j in 1:50)
            lines(object$procBeq0[,1],as.matrix(object$sim.test.procBeq0[[j]])[,i-1],
                  col="grey",lwd=1,lty=1,type="s")
          lines(object$procBeq0[,1],object$procBeq0[,i],lwd=2,type="s")
	}  ## }}}

    } ## }}} 
    else if (score==2)  ## {{{  plot score proces
    {
      dim1<-length(object$obs.test)
      if (sum(specific.comps)==FALSE) comp<-1:dim1 else comp<-specific.comps

      if (!is.null(xlab)) { if (length(xlab)==1) xlab <- rep(xlab,length(comp)); } 
      if (!is.null(main)) { if (length(main)==1) main <- rep(main,length(comp)); mains <- FALSE; } 

      v <- 0
      for (i in comp) ## {{{ 
      {
	v <- v+1
        if (nrow(object$proc.cumz[[i]])==1) TYPE<-"p" else TYPE<-"l"; 

        if (is.null(xlab)) xlabl <- colnames(object$proc.cumz[[v]])[1] else xlabl <- xlab[i]
###	print(colnames(object$proc.cumz[[v]])); print(xlabl)

          if (TYPE=="l") 
            {
              ranyl<-range(object$proc.cumz[[i]][,2]);
              for (j in 1:50) ranyl<-range(c(ranyl,(object$sim.test.proccumz[[i]])[,j])); 
              mr<-max(abs(ranyl));

              if (add.to.plot==FALSE)
	      {
	        if (!is.null(ylim)) 
                plot(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],type=TYPE,
                ylim=ylim,lwd=2,xlab=xlabl,ylab=ylab,...)
	        else 
                plot(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],type=TYPE,
                ylim=c(-mr,mr),lwd=2,xlab=xlabl,ylab=ylab,...)
	      }
              else
                lines(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],type="l")

              if (!is.null(main)) title(main=main[i]); 
              if (mains==TRUE) title(main=colnames(object$proc.cumz[[i]])[1]); 

              if (TYPE=="l") for (j in 1:50)
                lines(object$proc.cumz[[i]][,1],object$sim.test.proccumz[[i]][,j],
                      col="grey",lwd=1,lty=1,type="l")
              if (TYPE=="p") for (j in 1:50)
                points(object$proc.cumz[[i]][,1],object$sim.test.proccumz[[i]][,j],pch=".")
              lines(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],lwd=2); 
            } 
	    ## Prediction bandds ## {{{
	    if (conf.band==TRUE) {
	     col.alpha<-0.2
	     col.ci<-"darkblue"
	     lty.ci<-2
	      if (col.alpha==0) col.trans <- col.ci
	      else
	      col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
	      if (level!=0.05) c.alpha<-percen(object$sim.test[,i],1-level)
	      else c.alpha<-object$conf.band.cumz[i];
	      t<-object$proc.cumz[[i]][,1]
	      ci<-c.alpha*object$robvar.cumz[1:length(t),i]^.5
	      #print(t); print(ci)
	      lines(t,ci , lwd=1, col=col.ci, lty=lty.ci)
	      lines(t,-ci , lwd=1, col=col.ci, lty=lty.ci)
	      tt <- c(t, rev(t))
	      yy <- c(ci, rev(-ci))
	      polygon(tt,yy, col=col.trans, lty=0)      
	  } ## }}}
  } ## }}}

    } ## }}}

} ## }}}

"summary.cum.residuals" <- function (object,digits=3,...) 
{## {{{
  if (!inherits(object, 'cum.residuals')) stop ("Must be an cum.residuals object")

  # We print information about object:  
  cat("Test for cumulative MG-residuals \n\n")

  mtest<-(sum(object$conf.band)>0)

  if (mtest==FALSE) { 
    cat("Grouped cumulative residuals not computed, you must provide\n")
    cat("modelmatrix to get these (see help) \n\n")
  }
  if (mtest==TRUE) { 
    test0<-cbind(object$obs.testBeq0,object$pval.testBeq0)
    test0.is<-cbind(object$obs.testBeq0.is,object$pval.testBeq0.is) 
    colnames(test0)<- c("sup|  hat B(t) |","p-value H_0: B(t)=0")
    colnames(test0.is)<- c("int ( B(t) )^2 dt","p-value H_0: B(t)=0")  

    cat("Grouped Residuals consistent with model \n\n")
    prmatrix(round(test0,digits))
    cat("\n")
    prmatrix(round(test0.is,digits))
    cat("\n")
  }

  cumtest<-!is.null(object$obs.test)
  if (cumtest==FALSE) { 
    cat("Cumulative tests versus covariates not computed \n\n")
    cat("cum.resid=1 to compute these \n\n"); 
  }
  if (cumtest==TRUE) { 
    test0<-cbind(object$obs.test,object$pval.test)
    colnames(test0)<- c("sup|  hat B(t) |","p-value H_0: B(t)=0")

    cat("Residual versus covariates consistent with model \n\n")
    prmatrix(round(test0,digits))
  }
###  cat("  \n");cat("  Call: \n");dput(attr(object, "Call"));
  cat("\n");
} ## }}}


