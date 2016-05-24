aalen<-function (formula = formula(data),
     data = sys.parent(), start.time = 0, max.time = NULL, 
     robust=1, id=NULL, clusters=NULL, residuals = 0, n.sim = 1000,  
     weighted.test= 0,covariance=0,resample.iid=0,
     deltaweight=1,silent=1,weights=NULL,max.clust=1000,
     gamma=NULL,offsets=0){ ## {{{
## {{{ setting up variables 
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) { robust <- 1;}
  if (covariance==1 & robust==0) { covariance<-0;}
  if (sim==1 & robust==0) { n.sim<-0;sim <- 0}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
    m$start.time <- m$weighted.test <- m$max.time <- m$robust <- 
    m$weights <- m$residuals <- m$n.sim <- m$id <- m$covariance <- 
    m$resample.iid <- m$clusters <- m$deltaweight<-m$silent <- 
    m$max.clust <- m$gamma <- m$offsets <- NULL
  special <- c("const","cluster") 
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
  if(is.null(clusters)) clusters <- des$clusters ##########
  pxz <- px + pz; 

  if ((nrow(X)!=nrow(data) && (!is.null(id)))) stop("Missing values in design matrix not allowed with id \n"); 
###  if (nrow(Z)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 
  if (!is.null(id)) {
	  if (length(id)!=nrow(X)) stop("id length and data not the same\n"); 
  }

  cluster.call<- clusters; 
  survs<-read.surv(m,id,npar,clusters,start.time,max.time,silent=silent)
  times<-survs$times; 
  id<-survs$id.cal; 
  id.call<-id; 
  clusters<-gclusters <- survs$clusters; 
  stop.call <- time2<-survs$stop
  start.call <- survs$start
  status<-survs$status; 
  orig.max.clust <- survs$antclust
  dtimes <- sort(survs$stop[survs$status==1])
  nobs <- nrow(X); 
  if (is.null(weights)) weights <- rep(1,nrow(X)); 
  weights.call <- weights; 

  if ((!is.null(max.clust))) if (max.clust<survs$antclust) {
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- gclusters <-  as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
###	clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
	survs$antclust <- max.clust    
  }                                                         

  if ( (attr(m[, 1], "type") == "right" ) ) {  ## {{{
   ot <- order(-time2,status==1); 
   time2<-time2[ot]; 
   status<-status[ot]; 
   X<-as.matrix(X[ot,])
   if (npar==FALSE) Z<-as.matrix(Z[ot,])
   survs$stop<-time2;
   clusters<-clusters[ot]
   id<-id[ot];
   entry=rep(-1,nobs); 
   weights <- weights[ot]
   if (sum(offsets)!=0) offsets <- offsets[ot]
  } else {
        eventtms <- c(survs$start,time2)
        status <- c(rep(0, nobs), status)
        ix <- order(-eventtms,status==1)
        etimes    <- eventtms[ix]  # Entry/exit times
	status <- status[ix]
        survs$stop  <- etimes; 
        survs$start <- c(survs$start,survs$start)[ix]; 
        tdiff    <- c(-diff(etimes),start.time) # Event time differences
        entry  <- c(rep(c(1, -1), each = nobs))[ix]
	weights <- rep(weights, 2)[ix]
	X        <- X[rep(1:nobs, 2)[ix],]
	if (npar==FALSE) Z <- Z[rep(1:nobs,2)[ix],]
	id <- rep(id,2)[ix]
	clusters <- rep(clusters,2)[ix]
	if (sum(offsets)!=0) offsets <- rep(offsets,2)[ix]
} ## }}}

ldata<-list(start=survs$start,stop=survs$stop,antpers=survs$antpers,antclust=survs$antclust);

## }}}

  if (npar== TRUE) { ## {{{ Aalen model 
   ud <- aalenBase(times, ldata, X, status, id, clusters, robust = robust, 
   sim = sim, retur = residuals, antsim = n.sim,
   weighted.test = weighted.test,covariance=covariance,
   resample.iid=resample.iid,namesX=covnamesX,
   silent=silent,weights=weights,entry=entry,offsets=offsets)

    colnames(ud$cum) <- colnames(ud$var.cum) <- c("time",covnamesX)
    if (robust == 1) colnames(ud$robvar.cum) <- c("time", covnamesX)
    if (sim >= 1) {
      colnames(ud$test.procBeqC) <- c("time", covnamesX)
      names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- names(ud$obs.testBeqC.is) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC) <- colnames(ud$sim.testBeqC.is) <- covnamesX
      ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
    }
  } ## }}}
  else { ## {{{ Semiparametric additive risk model 
    if (px == 0) stop("No nonparametric terms (needs one!)")
    ud<-semiaalen(times, ldata, X, Z, 
    status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
    weighted.test = weighted.test, retur =
	residuals,covariance=covariance,
	resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
	deltaweight=deltaweight,gamma=gamma,
	silent=silent,weights=weights,entry=entry,offsets=offsets)

    if (px > 0) {
      colnames(ud$cum) <- colnames(ud$var.cum) <- c("time", covnamesX)
      if (robust == 1) 
        colnames(ud$robvar.cum) <- c("time", covnamesX)
      if (sim >= 1) {
        colnames(ud$test.procBeqC) <- c("time", covnamesX)
        names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC.is) <- colnames(ud$sim.testBeqC) <- covnamesX
        ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
      }
    }
    ud$gamma<-as.matrix(ud$gamma);
    rownames(ud$gamma) <- c(covnamesZ)
    rownames(ud$intZHdN) <- c(covnamesZ)
    colnames(ud$gamma) <- "estimate"
    colnames(ud$var.gamma) <- c(covnamesZ)
    rownames(ud$var.gamma) <- c(covnamesZ)
    colnames(ud$robvar.gamma) <- c(covnamesZ)
    colnames(ud$intZHZ) <- c(covnamesZ)
    rownames(ud$var.gamma) <- c(covnamesZ)
  } ## }}}

  attr(ud,"stratum")<-ud$stratum; 
  attr(ud, "Call") <- call
  attr(ud, "Formula") <- formula
  attr(ud, "id") <- id.call
  attr(ud, "cluster.call") <- cluster.call
  attr(ud, "cluster") <- gclusters
  attr(ud, "start.time") <- start.time
  attr(ud, "stop") <- stop.call
  attr(ud, "start") <- start.call
  attr(ud, "status") <- survs$status
  attr(ud, "residuals") <- residuals
  attr(ud, "max.clust") <- max.clust; 
  attr(ud, "max.time") <- max.time; 
  attr(ud, "weights") <- weights.call; 
  attr(ud, "orig.max.clust") <- orig.max.clust 
  class(ud) <- "aalen"
  ud$call<-call
  return(ud)
} ## }}}

plot.aalen <-  function (x, pointwise.ci=1, hw.ci=0,
sim.ci=0, robust=0, specific.comps=FALSE,level=0.05, start.time = 0, 
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE,...) 
{ ## {{{
  object <- x; rm(x);
  if (!inherits(object,'aalen') ) 
    stop ("Must be output from Aalen function") 
 
  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci, 
        hw.ci=hw.ci,
        sim.ci=sim.ci, robust=robust, specific.comps=specific.comps,level=level, 
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot, 
        mains=mains, xlab=xlab, ylab =ylab) 
  else plotScore(object, specific.comps=specific.comps, mains=mains,
                  xlab=xlab,ylab =ylab); 
} ## }}}

"print.aalen" <- function (x,...) 
{ ## {{{
summary.aalen(x,...)
###  object <- x; rm(x);
###  if (!inherits(object, 'aalen')) stop ("Must be an aalen object")
###
###  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
###    
###                                        # We print information about object:
###  
###  cat("Additive Aalen Model \n\n")
###  cat(" Nonparametric terms : "); cat(colnames(object$cum)[-1]); cat("   \n");  
###  if (semi) {
###    cat(" Parametric terms :  "); cat(rownames(object$gamma)); 
###    cat("   \n");  } 
###  cat("   \n");  
###
###  cat("  Call: \n"); dput(attr(object, "Call")); cat("\n"); 
} ## }}}

"summary.aalen" <-
function (object,digits = 3,...) 
{ ## {{{
  aalen.object <- object; rm(object);
  
  obj<-aalen.object
  if (!inherits(aalen.object, 'aalen')) stop ("Must be an aalen object")
  
  if (is.null(aalen.object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Additive Aalen Model \n\n")
 #cat("Nonparametric terms : "); cat(colnames(aalen.object$cum)[-1]);
 #cat("   \n");  

  timetest(obj,digits=digits); 

  if (semi) {
    cat("Parametric terms :  ");  #cat(rownames(aalen.object$gamma)); 
  }
  cat("   \n");  

  if (semi) {
    out=coef.aalen(aalen.object,digits=digits); 
    out=signif(out,digits=digits)
    print(out)

    if (is.null(aalen.object$pstest.pval)==FALSE) {
      res<-cbind(aalen.object$sup.pscore,aalen.object$pstest.pval);
      colnames(res)<-c("sup of pseudo-score test","p-value H_0: B(t)=b t");  
      cat(" \n");  
      cat("Test for time invariant effects \n")
      prmatrix(signif(res,digits))
    }
  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(aalen.object, "Call"))
  cat("\n")
} ## }}}

coef.aalen <- function(object, digits=3,...) {
   coefBase(object,digits=digits,...)
}
