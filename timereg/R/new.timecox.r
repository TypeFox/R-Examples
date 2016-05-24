const<-function(x) x

timecox<-function(formula=formula(data),data=sys.parent(),
start.time=0,max.time=NULL,id=NULL,clusters=NULL,
n.sim=1000,residuals=0,robust=1,Nit=20,bandwidth=0.5,
method="basic",weighted.test=0,degree=1,covariance=0)
{
  sim2<-0; if (n.sim==0) sim<-0 else sim<-1;
  if (method!="basic") stop("Only runs the default method at the moment\n"); 

                                        #if (resample.iid==1 & robust==0) {
                                        #cat("When robust=0 no iid representaion computed\n"); 
                                        #resample.iid<-0;}
  if (covariance==1 & robust==0) {
    cat("When robust=0 no covariance computed \n"); 
    cat("covariance set to 0\n"); 
    covariance<-0;}
  if (sim==1 & robust==0) {
    cat("When robust=0, No simulations \n"); 
    cat("n.sim set to 0\n"); 
    n.sim<-0;}
  if (residuals==1 & robust==0) {
    cat("When robust=0, no martingale residuals \n"); 
    cat("residuals set to 0\n"); 
    residuals<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}

  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$robust<-m$start.time<-m$degree<-m$weighted.test<-
    m$method<-m$Nit<-m$bandwidth<-m$max.time<-m$pers<-
      m$residuals<-m$n.sim<-m$id<-m$covariance<-NULL

  special <- c("const","cluster")
  Terms <- if(missing(data)) terms(formula, special)
  else              terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  
  
  if (is.null(Z)==TRUE) XZ<-X else XZ<-cbind(X,Z); 

  if (method=="breslow" && intercept==1) { 
    covnamesX<-covnamesX[-1]; X<-as.matrix(X[,-1]); XZ<-as.matrix(XZ[,-1]); 
    colnames(X)<-covnamesX; px<-px-1;}
  pxz <- px + pz;

  survs<-read.surv(m,id,npar,clusters,start.time,max.time,model="timecox")
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; 
  time2<-survs$stop; time<-survs$start
  status<-survs$status; Ntimes<-sum(status); 
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  times<-c(start.time,time2[status==1]); times<-sort(times);
  Ntimes <- Ntimes+1; 
  if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
  times<-times[times<maxtimes];
  bandwidth<-(maxtimes-start.time)*bandwidth; 

###  if (sum(time)>0) stop("Delayed entry data not allowed for this function \n"); 

  if (method=="breslow") 
    beta<-coxph(Surv(time,time2,status)~XZ)$coef
  else if (method=="basic" && intercept==1) 
    beta<-coxph(Surv(time,time2,status)~XZ[,-1])$coef
  else beta<-coxph(Surv(time,time2,status)~XZ)$coef;
  beta0<-c(0,0,beta)  
  if (method=="basic" && intercept==0)  beta0<-c(0,beta); 
  bhat<-matrix(beta0,length(times),length(beta0),byrow=TRUE); 
  timerange<-range(times); 
  bhat[,1]<-times; 
  if (method=="breslow" || intercept==1) {
    bhat[,2]<-sum(status)/sum(ldata$stop-ldata$start);
    if (method=="basic") bhat[,2]<-log(bhat[,2]); 
  }

  if (npar==TRUE) {
    #cat("Nonparametric Multiplicative Hazard Model"); cat("\n");
    ud<-timecoxBase(times,ldata,X,status,id,bhat,
                    sim=sim,antsim=n.sim,degree=degree,robust=robust,
                    band=bandwidth,it=Nit,method=method,retur=residuals,sim2=sim2,
                    weighted.test=weighted.test,covariance=covariance);

    if (method=="breslow") covnamesX<-c("Cumulative Baseline",covnamesX); 
    colnames(ud$cum)<-colnames(ud$var.cum)<-c("time",covnamesX)
    if (robust==1) colnames(ud$robvar.cum)<-c("time",covnamesX)

    if (sim==1) {
      #if (method=="breslow") covnamesX<-covnamesX[-1]; 
      colnames(ud$test.procBeqC)<- c("time",covnamesX)
      names(ud$conf.band)<-names(ud$pval.testBeq0)<- 
        names(ud$pval.testBeqC)<- names(ud$pval.testBeqC.is)<- 
          names(ud$obs.testBeqC.is)<- 
            names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<- covnamesX; 
      colnames(ud$sim.testBeq0)<- colnames(ud$sim.testBeqC)<- 
        colnames(ud$sim.testBeqC.is)<- covnamesX; 
      ud$sim.testBeqC.is<-ud$sim.testBeqC<-NULL;
      if (method=="breslow" && sim2==1) 
        names(ud$pval.testBeqC.is1)<-names(ud$pval.testBeqC.is2)<-
          names(ud$obs.testBeqC.is1)<-names(ud$obs.testBeqC.is2)<- covnamesX; 
    }
  }
  else {
                                        #cat("Semiparametric Multiplicative Risk Model"); cat("\n");
    if (px==0) { stop("No nonparametric terms (needs one!)"); }

    if (method=="breslow") {
                                        #print(c(px,pxz)); 
      gamma<-bhat[1,(px+3):(pxz+2)]; bhat<-bhat[,1:(px+2)]; } 
    else {
      gamma<-bhat[1,(px+2):(pxz+1)]; bhat<-bhat[,1:(px+1)] }
                                        #print(gamma); print(bhat)
                                        #print(X[1:5,]); print(Z[1:5,]); 

    ud<-semicox(times,ldata,X,Z,
                status,id,bhat,gamma=gamma,sim=sim,antsim=n.sim,
                band=bandwidth,it=Nit,method=method,retur=residuals,robust=robust,
                degree=degree,weighted.test=weighted.test,covariance=covariance)

    if (px>0) {
      if (method=="breslow") 
        colnames(ud$cum)<- colnames(ud$var.cum)<- 
          c("time","Cumulative Baseline",covnamesX)
      else 
        colnames(ud$cum)<- colnames(ud$var.cum)<- c("time",covnamesX)
      if (robust==1) { 
        if (method=="breslow") colnames(ud$robvar.cum)<- 
          c("time","Cumulative Baseline",covnamesX) else
        colnames(ud$robvar.cum)<-c("time",covnamesX);  }
      if (sim>=1) {
        if (method=="breslow") name<-
          c("time","Cumulative Baseline",covnamesX) else name<-c("time",covnamesX)

        colnames(ud$test.procBeqC)<- name;
        names(ud$conf.band)<- names(ud$pval.testBeq0)<- 
          names(ud$pval.testBeqC)<- 
            names(ud$pval.testBeqC.is)<- names(ud$obs.testBeqC.is)<- 
              names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<- 
                colnames(ud$sim.testBeq0)<- colnames(ud$sim.testBeqC.is)<- 
                  colnames(ud$sim.testBeqC)<- name[-1]; 
        ud$sim.testBeqC.is<-ud$sim.testBeqC<-NULL;
      }
    }

    rownames(ud$gamma)<-c(covnamesZ); 
    colnames(ud$gamma)<-"estimate"; 
    colnames(ud$var.gamma)<-c(covnamesZ); 
    rownames(ud$var.gamma)<-c(covnamesZ); 
    colnames(ud$robvar.gamma)<-c(covnamesZ); 
    rownames(ud$var.gamma)<-c(covnamesZ); 
  }

  ud$method<-method
  attr(ud,"Call")<-call; 
  class(ud)<-"timecox"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster")<-cluster.call;
  attr(ud,"start.time") <- start.time
  attr(ud,"start")<- time;
  attr(ud,"stop")<- time2;
  attr(ud,"status")<-status;
  attr(ud,"time2")<-time2;
  attr(ud,"residuals")<-residuals;
  attr(ud,"max.time")<-max.time;
  attr(ud,"stratum")<-0;
  ud$call<-call
  return(ud); 
}

"plot.timecox" <-  function (x,..., pointwise.ci=1,
hw.ci=0, sim.ci=0, robust=0, specific.comps=FALSE,level=0.05,
start.time = 0,
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE)
{
  object <- x; rm(x);
  if (!inherits(object,'timecox') ) stop ("Must be output
from Cox-Aalen function")

  if (score==FALSE) plot.cums(object,
        pointwise.ci=pointwise.ci,
        hw.ci=hw.ci,
        sim.ci=sim.ci, robust=robust,
        specific.comps=specific.comps,level=level,
        start.time = start.time, stop.time = stop.time,
        add.to.plot=add.to.plot,
        mains=mains, xlab=xlab, ylab =ylab)
  else plotScore(object, specific.comps=specific.comps,
                  mains=mains,
                  xlab=xlab,ylab =ylab);
}

"print.timecox"<-
function (x,...) 
{
  timecox.object <- x; rm(x);
  if (!inherits(timecox.object, 'timecox')) 
    stop ("Must be an timecox.object")

if (is.null(timecox.object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Multiplicative Hazard Model \n\n")
  cat(" Nonparametric terms : "); cat(colnames(timecox.object$cum)[-1]);
  cat("   \n");  
  if (semi) {
  cat(" Parametric terms :  "); cat(rownames(timecox.object$gamma)); 
  cat("   \n");  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(timecox.object, "Call"))
  cat("\n")
}


"summary.timecox" <-
function (object,..., digits = 3) 
{
  timecox.object <- object; rm(object);
  obj<-timecox.object
  if (!inherits(timecox.object, 'timecox')) 
    stop ("Must be an timecox.object")
  
  if (is.null(timecox.object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
                                        # We print information about object:  
  cat("Multiplicative Hazard Model \n\n")

  timetest(obj,digits=digits); 

  if (obj$method=="breslow" && (!semi) && (obj$obs.testBeqC.is1!=FALSE)) {
    testsupBL<-cbind(obj$obs.testBeqC.is1,obj$pval.testBeqC.is1)
    testssBL<-cbind(obj$obs.testBeqC.is2,obj$pval.testBeqC.is2)
    cat("Tests without baseline correction\n")
    cat("BL(t) = int_0^t lambda_0(t) b(t) dt, L(t) = int_0^t lambda_0(t) dt \n")
    colnames(testsupBL)<-c("sup| BL(t) - (t/tau)B(tau) L(t)|","p-value H_0: B(t)=b t")  
    colnames(testssBL)<-c("int (BL(t)-(t/tau)B(tau) L(t))^2dt","p-value H_0: B(t)=b t")
    prmatrix(signif(testsupBL,digits))
    prmatrix(signif(testssBL,digits)) }


  if (semi) {
    cat("Parametric terms :  "); #cat(rownames(timecox.object$gamma)); 
  }
  cat("   \n");  

  if (semi) { 
    out=coef.timecox(timecox.object,digits=digits); 
    out=signif(out,digits=digits)
    print(out)
  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(timecox.object, "Call"))
  cat("\n")
}

coef.timecox<- function(object,..., digits=3) {
   coefBase(object,digits=digits)
}

