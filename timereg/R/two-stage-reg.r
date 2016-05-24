Gprop<-function(x) x

two.stage<-function(margsurv,data=sys.parent(),
Nit=60,detail=0,start.time=0,max.time=NULL,id=NULL,clusters=NULL,
robust=1,theta=NULL,theta.des=NULL,var.link=0,step=0.5,notaylor=0,se.clusters=NULL)
{ ## {{{
## {{{ seting up design and variables
rate.sim <- 1; 
if (class(margsurv)!="coxph") {
 formula<-attr(margsurv,"Formula");
 beta.fixed <- attr(margsurv,"beta.fixed")
 if (is.null(beta.fixed)) beta.fixed <- 1; 
 ldata<-aalen.des(formula,data=data,model="cox.aalen");
 id <- attr(margsurv,"id"); 
 mclusters <- attr(margsurv,"cluster")
 mclustind <- attr(margsurv,"cluster")
 cluster.call <- attr(margsurv,"cluster.call")
 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
 time2 <- attr(margsurv,"stop"); 
 start <- attr(margsurv,"start")
 antpers<-nrow(X);
 if (beta.fixed==1) Z <- NULL; 
 if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
 if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
 px<-ncol(X);
 if (is.null(clusters) && is.null(mclusters)) 
	 stop("No cluster variabel specified in marginal or twostage call \n"); 
 if (is.null(clusters)) { clusters <- mclusters; cluster.call <- cluster.call} else {cluster.call <- clusters;}

 if (is.null(se.clusters)) secluster <- mclusters;
 antsecluster <- length(unique(secluster))
 if (is.numeric(secluster)) secluster <-  sindex.prodlim(unique(secluster),secluster)-1 else  {
      seclusters <- as.integer(factor(clusters, labels = 1:antsecluster))-1
 }

### print("two-stage"); print(head(cluster.call))
 if (is.null(cluster.call)) notaylor <- 1
 if (is.null(margsurv$gamma.iid)) notaylor <- 1

} else { ## coxph ## {{{ 
  notaylor <- 1
  antpers <- margsurv$n
  id <- 0:(antpers-1)
  mt <- model.frame(margsurv)
  Y <- model.extract(mt, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
   if (attr(Y, "type") == "right") {
        time2 <- Y[, "time"]; 
        status <- Y[, "status"]
	start <- rep(0,antpers);
	} else {
	 start <- Y[, 1]; time2 <- Y[, 2];status <- Y[, 3];
        }
   Z <- matrix(1,antpers,length(coef(margsurv)));

   if (is.null(clusters)) stop("must give clusters for coxph\n");
   cluster.call <- clusters 
   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
   px <- 1; pz <- ncol(Z); 
   beta.fixed <- 0
   semi <- 1
   start.time <- 0
   if (is.null(se.clusters)) secluster <- clusters;
   antsecluster <- length(unique(secluster))
   if (is.numeric(secluster)) secluster <-  sindex.prodlim(unique(secluster),secluster)-1 else  {
       clusters <- as.integer(factor(clusters, labels = 1:antsecluster))-1
   }
}  ## }}} 

 if (length(clusters)!=length(secluster)) stop("length of se.clusters not consistent with cluster length\n"); 

  if (anyNA(clusters)) stop("Missing values in cluster varaibles\n"); 
  out.clust <- cluster.index(clusters);  
  clusters <- out.clust$clusters
  maxclust <- out.clust$maxclust 
  antclust <- out.clust$antclust
  idiclust <- out.clust$idclust
  cluster.size <- out.clust$cluster.size

  if (sum(abs(start))>0) lefttrunk <- 1  else lefttrunk <- 0;  
  cumhazleft <- 0; 
  RR <-  rep(1,antpers); 

  update <- 1;
  if (update==0) { ## {{{
  if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { ### compute cum hazards in time point infty; 
	  nn <- nrow(margsurv$cum) 
	  cum <- Cpred(margsurv$cum,time2)[,-1]
	  if (npar==TRUE) cumhaz <- apply(cum*X,1,sum)
	  if (npar==FALSE) cumhaz <- apply(cum*X,1,sum)*exp( Z %*% margsurv$gamma)
	  if (lefttrunk==1) {
	     cum <- Cpred(margsurv$cum,start)[,-1]
	     cumhazleft <- apply(cum*X,1,sum)
	     if (npar==TRUE) cumhazleft <-  cumhazleft
	     if (npar==FALSE) cumhazleft <- cumhazleft * exp( Z %*% margsurv$gamma)
	  } 
  } else { residuals<-margsurv$residuals$dM; cumhaz<-status-residuals; }
  } ## }}}

  if (update==1) 
  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
     if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { 
         resi <- residualsTimereg(margsurv,data=data) 
         residuals <- resi$residuals; 
	 cumhaz <- resi$cumhaz; 
	 cumhazleft <- resi$cumhazleft; 
	 RR <- resi$RR
     } else { residuals <- margsurv$residuals$dM; 
              cumhaz <- status-residuals; 
	      if (class(margsurv)=="cox.aalen") RR  <- exp( Z %*% margsurv$gamma)
     }
  }
  else if (class(margsurv)=="coxph") {
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       cumhazleft <- rep(0,antpers)
       RR<- exp(margsurv$linear.predictors-sum(margsurv$means*coef(margsurv)))
        if ((lefttrunk==1)) { 
           baseout <- basehaz(margsurv,centered=FALSE); 
           cum <- cbind(baseout$time,baseout$hazard)
	   cum <- Cpred(cum,start)[,2]
	   cumhazleft <- cum * RR 
	}
  } ## }}}

###  print(head(cbind(residuals,cumhaz,RR,time2,status)))

  ratesim<-rate.sim; 
  inverse<-var.link
  pxz <- px + pz;
  times<-c(start.time,time2[status==1]); 
  times<-sort(times);
  if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
  times<-times[times<maxtimes]
  Ntimes <- sum(status[time2<maxtimes])+1; 


  Biid<-c(); gamma.iid <- 0; 
  if (is.null(margsurv$B.iid)) notaylor <- 1; 
  if (notaylor==0) {
    nBiid <- length(margsurv$B.iid)
    if (nBiid!=antsecluster) stop("Number of clusters for marginal models must be consistent with se.cluster for standard error sandwich\n"); 
    for (i in 1:nBiid) Biid<-cbind(Biid,margsurv$B.iid[[i]]); 
    if (!is.null(margsurv$gamma.iid)) gamma.iid<-margsurv$gamma.iid;
    if (is.null(margsurv$time.sim.resolution)) { 
	   time.group <- (1:nrow(Biid))-1; 
           maxtimesim <- nrow(Biid); 
	   timereso <- margsurv$cum[,1] 
    }  
    else {
      timereso <- margsurv$time.sim.resolution
       qqc <- cut(times, breaks = margsurv$time.sim.resolution, include.lowest = TRUE)    
       time.group <- as.integer(factor(qqc, labels = 1:(nrow(Biid)-1)))
       maxtimesim <- nrow(Biid); 
    } 
    if (class(margsurv)=="cox.aalen")  { 
       times <- margsurv$time.sim.resolution 
       Ntimes <-length(times)
    }
  } else {time.group <- 1; maxtimesim <- 1; timereso <- 1}

  if (is.null(theta.des)==TRUE) ptheta<-1; 
  if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else
  theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des); 
  if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE) {
      if (var.link==1) theta<- rep(-1,ptheta); 
      if (var.link==0) theta<- rep(exp(-5),ptheta); 
  }
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters !, maxclust size=1\n"); 
  theta.iid <- matrix(0,antsecluster,ptheta)
  ## }}}

  
  nparout <- .C("twostagereg", 
        as.double(times), as.integer(Ntimes), as.double(X),
   	as.integer(antpers), as.integer(px), as.double(Z), 
	as.integer(antpers), as.integer(pz), as.integer(antpers),         ## 9 
	as.double(start),as.double(time2), as.integer(Nit), 
	as.integer(detail), as.integer(id), as.integer(status),           ## 15
	as.integer(ratesim), as.integer(robust), as.integer(clusters),    
	as.integer(antclust), as.integer(beta.fixed), as.double(theta),
	as.double(var.theta), as.double(theta.score), as.integer(inverse), 
	as.integer(cluster.size),                                          ## 25
	as.double(theta.des), as.integer(ptheta), as.double(Stheta),
	as.double(step), as.integer(idiclust), as.integer(notaylor),
	as.double(gamma.iid),as.double(Biid),as.integer(semi), as.double(cumhaz) ,
	as.double(cumhazleft),as.integer(lefttrunk),as.double(RR),
	as.integer(maxtimesim),as.integer(time.group),as.integer(secluster),
	as.integer(antsecluster),as.double(theta.iid), as.double(timereso),PACKAGE = "timereg")

## {{{ handling output
   gamma <- margsurv$gamma
   Varbeta <- margsurv$var.gamma; RVarbeta <- margsurv$robvar.gamma;
   score <- margsurv$score; Iinv <- margsurv$D2linv;
   cumint <- margsurv$cum; vcum <- margsurv$var.cum; Rvcu <- margsurv$robvar.cum;

   theta<-matrix(nparout[[21]],ptheta,1);  
   var.theta<-matrix(nparout[[22]],ptheta,ptheta); 
   theta.score<-nparout[[23]]; 
   SthetaI<-matrix(nparout[[28]],ptheta,ptheta); 
   theta.iid  <- matrix(nparout[[43]],antsecluster,ptheta); 
   theta.iid <- theta.iid %*% SthetaI

###  if (is.null(call.secluster) & is.null(max.clust)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)

   ud <- list(cum = cumint, var.cum = vcum, robvar.cum = Rvcu, 
       gamma = gamma, var.gamma = Varbeta, robvar.gamma = RVarbeta, 
       D2linv = Iinv, score = score,  theta=theta,var.theta=var.theta,
       SthetaInv=SthetaI,theta.score=theta.score,theta.iid=theta.iid)

  ptheta<-length(ud$theta); 
  if (ptheta>1) {
                rownames(ud$theta)<-colnames(theta.des);
                names(ud$theta.score)<-colnames(theta.des); } else { 
		names(ud$theta.score)<- rownames(ud$theta)<-"intercept" } 


  attr(ud,"Call")<-call; 
  class(ud)<-"two.stage"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id;
  attr(ud,"cluster")<-clusters;
  attr(ud,"cluster.call")<-cluster.call;
  attr(ud,"secluster")<-secluster;
  attr(ud,"start")<-start; 
  attr(ud,"time2")<-time2; 
  attr(ud,"var.link")<-var.link
  attr(ud,"beta.fixed")<-beta.fixed
  attr(ud,"marg.model")<-class(margsurv)

  return(ud) 
  ## }}} 
} ## }}} 
  
  summary.two.stage<-function (object,digits = 3,...) { ## {{{ 

  if (!(inherits(object, 'two.stage') )) stop("Must be a Two-Stage object")
  prop<-TRUE; 
  if (is.null(object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE
    
  var.link<-attr(object,"var.link");
  cat("Dependence parameter for Clayton-Oakes-Glidden  model\n"); 

  if (sum(abs(object$theta.score)>0.000001) ) 
    cat("Variance parameters did not converge, allow more iterations\n\n"); 

  ptheta<-nrow(object$theta)
  sdtheta<-diag(object$var.theta)^.5
  if (var.link==0) {
      vari<-object$theta
      sdvar<-diag(object$var.theta)^.5
  }
  else {
      vari<-exp(object$theta)
      sdvar<-vari*diag(object$var.theta)^.5
  }
  dep<-cbind(object$theta[,1],sdtheta)
  walddep<-object$theta[,1]/sdtheta; 
  waldpdep<-(1-pnorm(abs(walddep)))*2

  kendall<-1/(1+2/vari) 
  kendall.ll<-1/(1+2/(object$theta+1.96*sdvar)) 
  kendall.ul<-1/(1+2/(object$theta-1.96*sdvar)) 
  if (var.link==0) resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,kendall)),digits)
  else resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,vari,sdvar,kendall)),digits);

  if (var.link==0) colnames(resdep) <- c("Variance","SE","z","P-val","Kendall's tau") 
  else colnames(resdep)<-c("log(Variance)","SE","z","P-val","Variance","SE Var.",
                           "Kendall's tau")
  prmatrix(resdep); cat("   \n");  

  if (attr(object,"marg.model")!="coxph")
  if (attr(object,"beta.fixed")==0) { ## {{{ 
###  cat("Marginal Cox-Aalen model fit\n\n"); 
  if (sum(abs(object$score)>0.000001) && sum(object$gamma)!=0) 
    cat("Marginal model did not converge, allow more iterations\n\n"); 
###  if (prop) {
###    if (p.o==FALSE) cat("Proportional Cox terms :  \n") else  cat("Covariate effects \n")
###
###    out=coef.two.stage(object,digits=digits);
###    out=signif(out,digits=digits)
###    print(out)
###
###  }

  } ## }}} 
###   cat("   \n");  cat("  Call: \n"); dput(attr(object, "Call")); 
  cat("\n");
} ## }}}

print.two.stage <- function (x,digits = 3,...) { ## {{{
	summary.two.stage(x,digits=digits,...)
###  if (!(inherits(x, 'two.stage') )) stop("Must be a Two-Stage object")
###  cat(" Two-stage estimation for Clayton-Oakes-Glidden  model\n"); 
###  cat(" Marginals of Cox-Aalen form, dependence by variance of Gamma distribution\n\n");  
###  object <- x; rm(x);
###  
###  cat(" Nonparametric components : "); 
###  cat(colnames(object$cum)[-1]); cat("   \n");  
###  if (!is.null(object$gamma)) {
###    cat(" Parametric components :  "); cat(rownames(object$gamma)); 
###    cat("   \n");
###  } 
###  cat("   \n");  
###
###  cat(" Call: \n");
###  print(attr(object,'Call'))
} ## }}}

coef.two.stage<-function(object,digits=3,d2logl=1,...) { ## {{{ 

  if (!(inherits(object, 'two.stage') )) stop("Must be a Two-Stage object")
  var.link <- attr(object,"var.link")
  ptheta<-nrow(object$theta)
  sdtheta<-diag(object$var.theta)^.5

  if (var.link==0) {
      vari<-object$theta
      sdvar<-diag(object$var.theta)^.5
  }
  else {
      vari<-exp(object$theta)
      sdvar<-vari*diag(object$var.theta)^.5
  }
  dep<-cbind(object$theta[,1],sdtheta)
  walddep<-object$theta[,1]/sdtheta; 
  waldpdep<-(1-pnorm(abs(walddep)))*2

  kendall<-1/(1+2/vari) 
  kendall.ll<-1/(1+2/(object$theta+1.96*sdvar)) 
  kendall.ul<-1/(1+2/(object$theta-1.96*sdvar)) 
  if (var.link==0) resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,kendall)),digits)
  else resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,vari,sdvar,kendall)),digits);

  if (var.link==0) colnames(resdep) <- c("Variance","SE","z","P-val","Kendall's tau") 
  else colnames(resdep)<-c("log(Variance)","SE","z","P-val","Variance","SE Var.","Kendall's tau")
###  prmatrix(resdep); cat("   \n");  
  return(resdep)
} ## }}} 

plot.two.stage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
		level=0.05, 
		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
                xlab="Time",ylab ="Cumulative regression function",...) 
{ ## {{{
  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
  object <- x; rm(x);  
 
  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
  if (robust>=1) V<-object$robvar.cum; 

  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
  if (stop.time==0) stop.time<-max(B[,1]);

  med<-B[,1]<=stop.time & B[,1]>=start.time
  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
  Vrob<-object$robvar.cum; 
  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

  c.alpha<- qnorm(1-level/2)
  for (v in comp) { 
    c.alpha<- qnorm(1-level/2)
    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
    if (add.to.plot==FALSE) 
      {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
        if (mains==TRUE) title(main=colnames(B)[v]); }
    else lines(B[,1],est,type="s"); 
    if (pointwise.ci>=1) {
      lines(B[,1],ul,lty=pointwise.ci,type="s");
      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
    if (robust>=1) {
      lines(B[,1],ul,lty=robust,type="s"); 
      lines(B[,1],nl,lty=robust,type="s"); }
    abline(h=0); 
  }
}  ## }}}

predict.two.stage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,X2=NULL,Z2=NULL,
			      theta=NULL,theta.des=NULL,diag=TRUE,...)
{ ## {{{
time.coef <- data.frame(object$cum)
if (!is.null(times))   cum <- Cpred(object$cum,times)  else cum <- object$cum; 
if (!is.null(times2)) cum2 <- Cpred(object$cum,times2) else cum2 <- object$cum;

if (is.null(Z) & (!is.null(object$gamma))) Z <- matrix(0,1,nrow(object$gamma));
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}
if (is.null(X)) X <- 1;

if (is.null(X2)) X2 <- X;
if (is.null(Z2)) Z2 <- Z2;
if (is.null(Z2) & (!is.null(object$gamma))) Z2 <- Z
if (is.null(X2) & (!is.null(Z2)))  {Z2 <- as.matrix(Z2);  X2 <- matrix(1,nrow(Z2),1)}
if (is.null(Z2) & (!is.null(X2)))  {X2 <- as.matrix(X2);  Z2 <- matrix(0,nrow(X2),1); gamma <- 0}

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1]) 
   time.part2 <-  X2 %*% t(cum2[,-1]) 
   if (!is.null(object$gamma)) { 
	   gamma <- object$gamma
	   RR <- exp( Z %*% gamma ); 
	   RR2 <- exp( Z2 %*% gamma ); 
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR2 )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else { 
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
	time.part2 <-  apply(as.matrix(X2*cum2[,-1]),1,sum) 
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%object$gamma); 
	RR2 <- exp(Z2 %*%object$gamma); 
	cumhaz <- c((time.part) * RR) ;  
	cumhaz2 <- c((time.part2) * RR2); 
} else {
	cumhaz <- c(time.part);  cumhaz2 <- c(time.part2); 
} 
S1 <- pmin(1,exp(-cumhaz)); S2 <- pmin(1,exp(-cumhaz2))
###print(length(S1))
###print(length(S2))

if (is.null(theta))  theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% theta)
if (attr(object,"var.link")==1) theta  <- exp(theta) 

if (diag==FALSE) St1t2<- (outer(c(S1)^{-(theta)},c(S2)^{-(theta)},FUN="+") - 1)^(-(1/theta)) else 
St1t2<- ((S1^{-(theta)}+S2^{-(theta)})-1)^(-(1/theta))
###St1t2<- ((S1^{-(1/theta)}+S2^{-(1/theta)})-1)^(-(theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} ## }}}

