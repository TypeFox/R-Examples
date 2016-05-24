##' Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data 
##' using marginals that are on Cox or addtive form. 
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##'
##' @export
##' @references
##' Clayton-Oakes and Plackett bivariate survival distributions,
##' @examples
##' data(diabetes)
##' 
##' # Marginal Cox model  with treat as covariate
##' margph <- coxph(Surv(time,status)~treat,data=diabetes)
##' ### Clayton-Oakes, from timereg
##' fitco1<-two.stage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id)
##' summary(fitco1)
##' ### Plackett model
##' fitp<-twostage(margph,data=diabetes,theta=3.0,Nit=40,
##'                clusters=diabetes$id,var.link=1)
##' summary(fitp)
##' ### Clayton-Oakes
##' fitco2<-twostage(margph,data=diabetes,theta=0.0,detail=0,
##'                  clusters=diabetes$id,var.link=1,model="clayton.oakes")
##' summary(fitco2)
##' fitco3<-twostage(margph,data=diabetes,theta=1.0,detail=0,
##'                  clusters=diabetes$id,var.link=0,model="clayton.oakes")
##' summary(fitco3)
##'
##' ### without covariates using Aalen for marginals
##' marg <- aalen(Surv(time,status)~+1,data=diabetes,n.sim=0,max.clust=NULL,robust=0)
##' fitpa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,
##'                 clusters=diabetes$id,score.method="optimize")
##' summary(fitpa)
##' 
##' fitcoa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,
##'                  var.link=1,model="clayton.oakes")
##' summary(fitcoa)
##' 
##' ### Piecewise constant cross hazards ratio modelling
##' ########################################################
##' 
##' d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
##' udp <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",
##'                           id="cluster",timevar="time",
##'                           status="status",model="clayton.oakes",silent=0)
##' summary(udp)
##' 
##' \donttest{ ## Reduce Ex.Timings
##' ### Same model using the strata option, a bit slower
##' ########################################################
##' ## makes the survival pieces for different areas in the plane 
##' ##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##' ##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")
##' 
##' ## everything done in one call 
##' ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
##' ud$strata <- factor(ud$strata); 
##' ud$intstrata <- factor(ud$intstrata)
##' 
##' ## makes strata specific id variable to identify pairs within strata
##' ## se's computed based on the id variable across strata "cluster"
##' ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000
##' 
##' marg2 <- aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
##'                data=ud,n.sim=0,robust=0)
##' tdes <- model.matrix(~-1+factor(strata),data=ud)
##' fitp2<-twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes,step=0.5)
##' summary(fitp2)
##' 
##' ### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
##' ud$stratas <- ud$strata; 
##' ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
##' tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
##' fitp3<-twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes2,step=0.5)
##' summary(fitp3)
##' 
##' ### same model using strata option, a bit slower 
##' fitp4<-twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes2,step=0.5,strata=ud$strata)
##' summary(fitp4)
##' }
##' @keywords survival
##' @author Thomas Scheike
##' @export
##' @param margsurv Marginal model 
##' @param data data frame
##' @param score.method Scoring method
##' @param Nit Number of iterations
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.des Variance component design
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param notaylor Taylor expansion
##' @param model model
##' @param marginal.trunc marginal left truncation probabilities
##' @param marginal.survival optional vector of marginal survival probabilities 
##' @param marginal.status related to marginal survival probabilities 
##' @param strata strata for fitting, see example
##' @param se.clusters for clusters for se calculation with iid
##' @param max.clust max se.clusters for se calculation with iid
##' @param numDeriv to get numDeriv version of second derivative, otherwise uses sum of squared score 
twostage <- function(margsurv,data=sys.parent(),score.method="nlminb",Nit=60,detail=0,clusters=NULL,
		     silent=1,weights=NULL, control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=1,
                     step=0.5,notaylor=0,model="plackett",
		     marginal.trunc=NULL,
		     marginal.survival=NULL,
		     marginal.status=NULL,strata=NULL,
		     se.clusters=NULL,max.clust=NULL,numDeriv=1)
{ ## {{{
## {{{ seting up design and variables
rate.sim <- 1; sym=1; 
if (model=="clayton.oakes") dep.model <- 1
else if (model=="plackett") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n"); 
start.time <- NULL

if (!is.null(margsurv)) 
if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen") { ## {{{
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1; 
	 ldata<-aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id"); 
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
	 time2 <- attr(margsurv,"stop"); start.time <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n"); 
         if (is.null(clusters)) clusters <- mclusters else if (sum(abs(clusters-mclusters))>0) 
         cat("Warning: Clusters for marginal model different than those specified for two.stage\n"); 
         if (!is.null(attr(margsurv,"max.clust")))
         if ((attr(margsurv,"max.clust")< attr(margsurv,"orig.max.clust")) && (!is.null(mclusters))) 
		  cat("Warning: Probably want to estimate marginal model with max.clust=NULL\n"); 
	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
## }}}
} else { ### coxph ## {{{
	  notaylor <- 1
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"]; 
	      status <- Y[, "status"]
		start.time <- rep(0,antpers);
		} else {
		 start.time <- Y[, 1]; 
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
###	   Z <- na.omit(model.matrix(margsurv)[,-1]) ## Discard intercept
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   cluster.call <- clusters; 
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z); 
	   if (sum(abs(start.time))>0) lefttrunk <- 1 else lefttrunk <- 0
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
###	   start.time <- 0
} ## }}}

  if (!is.null(start.time)) {
  if (sum(abs(start.time))>0) lefttrunk <- 1  else lefttrunk <- 0;  
  } else lefttrunk <- 0

if (!is.null(margsurv)) 
  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
         resi <- residualsTimereg(margsurv,data=data) 
         RR  <- resi$RR
	 psurvmarg <- exp(-resi$cumhaz); 
         ptrunc <- rep(1,length(psurvmarg)); 
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
  } ## }}}
  else if (class(margsurv)=="coxph") {  ## {{{
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz); 
       cumhazleft <- rep(0,antpers)
       ptrunc <- rep(1,length(psurvmarg)); 
       RR<- exp(margsurv$linear.predictors-sum(margsurv$means*coef(margsurv)))
        if ((lefttrunk==1)) { 
         baseout <- survival::basehaz(margsurv,centered=FALSE); 
         cum <- cbind(baseout$time,baseout$hazard)
	 cum <- Cpred(cum,start.time)[,2]
	 ptrunc <- exp(-cum * RR)
	}
  } ## }}}


  if (!is.null(marginal.survival) ) {
      if (!is.null(margsurv)) {  
	  cat("must give either margsurv model or marginal.survival=probabilities,\n"); 
          cat(" uses marginal.survival\n"); 
      }
	if (lefttrunk==1)  cat("Warnings specify only your own survival weights for right-censored data\n"); 
###        if (length(marginal.survival)!=length(start.time)) stop(paste("marginal.survival must have length=",antpers,"\n"));  
        psurvmarg <- marginal.survival
	antpers <- length(marginal.survival)
        RR <-  rep(1,antpers); 
	if (!is.null(marginal.trunc)) ptrunc <- marginal.trunc else ptrunc <- rep(1,antpers);
	if (!is.null(marginal.status)) status <- marginal.status else stop("must give censoring status\n"); 
  } 

  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
  if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
  if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

  cluster.call <- clusters

  out.clust <- cluster.index(clusters);  
  clusters <- out.clust$clusters
  maxclust <- out.clust$maxclust 
  antclust <- out.clust$antclust
  clusterindex <- out.clust$idclust
  clustsize <- out.clust$cluster.size
  call.secluster <- se.clusters
  if (is.null(se.clusters)) { se.clusters <- clusters; antiid <- nrow(clusterindex);} else  {
      iids <-  unique(se.clusters); 
      antiid <- length(iids); 
      if (is.numeric(se.clusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
       else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
  }
  if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n"); 

  if ((!is.null(max.clust))) if (max.clust< antiid) {
        coarse.clust <- TRUE
	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
	se.clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(se.clusters))
	maxclust <- max.clust    
	antiid <- max.clusters
  }                                                        

  ratesim<-rate.sim; 
###  pxz <- px + pz;
  if (is.null(theta.des)==TRUE) ptheta<-1; 
  if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else
  theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des); 
  if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE) {
         if (var.link==1) theta<- rep(-0.7,ptheta);  
         if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
  }       
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters, maxclust size=1\n"); 
  ## }}}

  loglike <- function(par) 
  { ## {{{
       Xtheta <- theta.des %*% matrix(c(par),ptheta,1); 
       DXtheta <- array(0,c(1,1,1));

###   dyn.load("twostage.so")

      outl<-.Call("twostageloglike", ## {{{
      icause=status,ipmargsurv=psurvmarg, 
      itheta=c(par),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
      ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
      itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid) 
      ## }}}

    if (detail==3) print(c(par,outl$loglike))

    attr(outl,"gradient") <-outl$score 
    if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==2) ret <- outl else ret <- outl$score
    return(ret)
  } ## }}}

  if (score.method=="optimize" && ptheta!=1) {cat("optimize only works for d==1, score.mehod set to nlminb \n"); score.method <- "nlminb";}

  theta.iid <- NULL
  logl <- NULL
  p <- theta
  if (score.method=="fisher.scoring") { ## {{{
    oout <- 2;  ### output control for obj
    if (Nit>0) 
    for (i in 1:Nit)
    {
        out <- loglike(p)
	hess <- -1* out$Dscore
	if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- hess 
        if (detail==1) {## {{{
          cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
          cat("theta:");print(c(p))
          cat("loglike:");cat(c(out$loglike),"\n"); 
          cat("score:");cat(c(out$score),"\n"); 
	  cat("hess:\n"); cat(out$Dscore,"\n"); 
        }## }}}
        delta <- hessi %*% out$score *step 
        p <- p+delta* step
        theta <- p; 
	if (is.nan(sum(out$score))) break; 
        if (sum(abs(out$score))<0.00001) break; 
        if (max(theta)>20) break; 
    }
    if (!is.nan(sum(p))) { 
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    out <- loglike(p) 
    logl <- out$loglike
    score1 <- score <- out$score
    hess1 <- hess <- -1* out$Dscore 
    if (iid==1) theta.iid <- out$theta.iid
    }
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3;  ## to get jacobian
      hess <- numDeriv::jacobian(loglike,p)
    if (detail==1 ) cat("numDeriv hessian stop\n"); 
    }
    if (detail==1 & Nit==0) {## {{{
          cat(paste("Fisher-Scoring ===================: final","\n")); 
          cat("theta:");print(c(p))
          cat("loglike:");cat(c(out$loglike),"\n"); 
          cat("score:");cat(c(out$score),"\n"); 
	  cat("hess:\n"); cat(out$Dscore,"\n"); 
    }## }}}
    if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- diag(nrow(hess))
    ## }}}
  } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
    oout <- 0; 
    tryCatch(opt <- nlminb(theta,loglike,control=control),error=function(x) NA)
    if (detail==1) print(opt); 
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    theta <- opt$par
    out <- loglike(opt$par)
    logl <- out$loglike
    score1 <- score <- out$score
    hess1 <- hess <- -1* out$Dscore
    if (iid==1) theta.iid <- out$theta.iid
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3; ## returns score 
      hess <- numDeriv::jacobian(loglike,opt$par)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
    }
    hessi <- lava::Inverse(hess); 
  ## }}}
  } else if (score.method=="optimize" && ptheta==1) { ## {{{  optimizer
    oout <- 0; 
    if (var.link==1) {mino <- -20; maxo <- 10;} else {mino <- 0.001; maxo <- 100;}
    tryCatch(opt <- optimize(loglike,c(mino,maxo)));
    if (detail==1) print(opt); 
    opt$par <- opt$minimum
    theta <- opt$par
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$par)
    logl <- out$loglike
    score1 <- score <- out$score
    hess1 <- hess <- -1* out$Dscore
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3;  ## to get jacobian
      hess <- numDeriv::jacobian(loglike,theta)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
    }
    hessi <- lava::Inverse(hess); 
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$estimate)
    logl <- out$loglike
    score1 <- out$score
    hess1 <- out$Dscore
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 

## {{{ handling output
  loglikeiid <- NULL
  robvar.theta <- NULL
  var.theta <- -1*hessi
  if (iid==1) {
     theta.iid <- out$theta.iid %*% hessi
     if (is.null(call.secluster) & is.null(max.clust)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)
     robvar.theta  <- (t(theta.iid) %*% theta.iid) 
     loglikeiid <- out$loglikeiid
  }
  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- rep("intercept",ptheta)
###  if (length(thetanames)==nrow(theta)) rownames(theta) <- thetanames
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,model=model,robvar.theta=robvar.theta,
             theta.iid=theta.iid,loglikeiid=loglikeiid,
	     thetanames=thetanames,loglike=-logl,score1=score1,Dscore=out$Dscore,margsurv=psurvmarg); 
  class(ud)<-"twostage" 
  attr(ud, "Formula") <- formula
  attr(ud, "clusters") <- clusters
  attr(ud, "cluster.call") <- cluster.call
  attr(ud, "secluster") <- c(se.clusters)
  attr(ud,"sym")<-sym; 
  attr(ud,"var.link")<-var.link; 
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  attr(ud, "Type") <- model
  attr(ud, "response") <- "survival"
  return(ud);
  ## }}}

} ## }}}

##' @export
summary.twostage <-function (object,digits = 3,silent=0,...) { ## {{{
  if (!(inherits(object,"twostage"))) stop("Must be a Two-Stage object")
  
  var.link<-attr(object,"var.link");
  if ((object$model=="plackett") & (silent==0)) cat("Dependence parameter for Plackett model \n"); 
  if (attr(object,"response")=="binomial") response <- "binomial" else response <- "survival"
  if ((object$model=="clayton.oakes") & (silent==0)) cat("Dependence parameter for Clayton-Oakes model \n"); 

  if ((sum(abs(object$score))>0.0001) & (silent==0))  {
	  cat("    Variance parameters did not converge, allow more iterations.\n"); 
	  cat(paste("    Score:",object$score,"  \n")); 
  }

  coefs <- coef.twostage(object,response=response,...);

  res <- list(estimates=coefs, type=attr(object,"Type"))
  class(res) <- "summary.twostage"
  res
} ## }}}

##' @export
coef.twostage <- function(object,var.link=NULL,response="survival",...)
{ ## {{{
  theta <- object$theta
  if (is.null(var.link))
     if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0
     else vlink <- var.link
  se<-diag(object$var.theta)^0.5
  res <- cbind(theta, se )
  wald <- theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  if (response=="survival") { 
       if (object$model=="plackett") {
       spearman <- alpha2spear(theta,link=vlink)
       Dspear <- numDeriv::jacobian(alpha2spear,theta,link=vlink) 
       var.spearman <- Dspear %*% object$var.theta %*%  Dspear
       se.spearman <- diag(var.spearman)^.5
       res <- as.matrix(cbind(res, wald, waldp,spearman,se.spearman))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  else colnames(res) <- c("Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
       }
       if (object$model=="clayton.oakes") {
       kendall <- alpha2kendall(theta,link=vlink)
       Dken <- numDeriv::jacobian(alpha2kendall,theta,link=vlink) 
       var.kendall<- Dken %*% object$var.theta %*%  Dken
       se.kendall <- diag(var.kendall)^.5
       res <- as.matrix(cbind(res, wald, waldp,kendall,se.kendall))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Kendall tau","SE")
       else colnames(res) <- c("Coef.", "SE","z", "P-val","Kendall tau","SE")
       if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
       }
  }
  return(res)
} ## }}}

##' @export
alpha2spear <- function(theta,link=1) { ## {{{ 
   if (link==1) theta <- exp(theta)
if (length(theta)>1) {
   out <- c()
   for (thet in theta) {
   if (thet!=1) out <- c(out,( (thet+1)/(thet-1) -2* thet* log(thet)/ (thet-1)^2))
   else out <- c(out,0)
   }
} else { if (theta!=1) out <- ( (theta+1)/(theta-1) -2* theta* log(theta)/ (theta-1)^2) }

return(out)
} ## }}} 

##' @export
alpha2kendall <- function(theta,link=0) {  ## {{{ 
   if (link==1) theta <- exp(theta)
   return(1/(1+2/theta)) 
} ## }}} 

##' @export
print.twostage<-function(x,digits=3,...)
{ ## {{{
  print(x$call); 
  cat("\n")
  print(summary(x)); 
} ## }}}

##' @export
plot.twostage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
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

##' @export
predict.twostage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,theta.des=NULL,diag=TRUE,...)
{ ## {{{
time.coef <- data.frame(object$cum)
if (!is.null(times)) {
cum <- Cpred(object$cum,times);
cum2 <- Cpred(object$cum,times);
} else { cum <- object$cum; cum2 <- object$cum }
if (!is.null(times2)) cum2 <- Cpred(object$cum,times2);

if (is.null(X)) X <- 1;
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1]) 
   time.part2 <-  X %*% t(cum2[,-1]) 
   if (!is.null(object$gamma)) { RR <- exp( Z %*% gamma ); 
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else { 
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
	time.part2 <-  apply(as.matrix(X*cum2[,-1]),1,sum) 
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%gamma); 
	cumhaz <- t( t(time.part) * RR );  
	cumhaz2 <- t( t(time.part2) * RR )} else {
		cumhaz <- time.part;  cumhaz2 <- time.part2; 
} 
S1 <- exp(-cumhaz); S2 <- exp(-cumhaz2)

if (attr(object,"var.link")==1) theta  <- exp(object$theta) else theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% object$theta)

if (diag==FALSE) St1t2<- (outer(c(S1)^{-(theta)},c(S2)^{-(theta)},FUN="+") - 1)^(-(1/theta)) else 
St1t2<- ((S1^{-(theta)}+S2^{-(theta)})-1)^(-(1/theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} ## }}}

##' @export
piecewise.twostage <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,
            score.method="optimize",Nit=100,detail=0,silent=1,weights=NULL,
            control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=1,step=0.5,model="plackett",data.return=0)
{ ## {{{
ud <- list()
if (missing(cut2)) cut2 <- cut1; 
nc1 <- length(cut1); nc2 <- length(cut2)
names1 <- names2 <- c()
theta.mat <- se.theta.mat <- cor.mat <- score.mat <- se.cor.mat <- matrix(0,nc1-1,nc2-1); 
clusters <- data[,id]
cluster.call <- clusters
idi <- unique(data[,id]); 
###print(head(idi))

## {{{ 
###   se.clusters=NULL,max.clust=1000,
###  evt saette cluster se max.clust paa 
###  if (is.null(se.clusters)) { se.clusters <- clusters; antiid <- nrow(clusterindex);} else  {
###      iids <-  unique(seclusters); 
###      antiid <- length(iids); 
###      if (is.numeric(seclusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
###       else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
###  }
###  if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n"); 
###
###  if ((!is.null(max.clust))) if (max.clust< antiid) {
###        coarse.clust <- TRUE
###	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
###	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
###	se.clusters <- as.integer(qqc)-1
###	max.clusters <- length(unique(se.clusters))
###	maxclust <- max.clust    
###	antiid <- max.clusters
###  }                                                         
## }}} 

if (iid==1) { theta.iid <- matrix(0,length(idi),(nc1-1)*(nc2-1));
              rownames(theta.iid) <- idi
            } else theta.iid <- NULL

k <- 0; 
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n");
datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
			status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent) 
if (silent<=-1) print("back in piecewise.twostage"); 
if (silent<=-1) print(summary(datalr)); 
if (silent<=-1) print(head(datalr)); 
if (silent<=-1) print(summary(datalr[,id])); 
 boxlr <- list(left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]))
### marg1 <- aalen(Surv(datalr$left,datalr[,timevar],datalr[,status])~+1,data=datalr,n.sim=0,max.clust=NULL,robust=0)
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
###

if (is.null(covars)) 
f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,")")))
else f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,"):",covars)))
marg1 <- aalen(f,data=datalr,n.sim=0,robust=0)

fitlr<-  twostage(marg1,data=datalr,clusters=datalr$tsid,model=model,score.method=score.method,
             Nit=Nit,detail=detail,silent=silent,weights=weights,
             control=control,theta=theta,theta.des=theta.des,var.link=var.link,iid=iid,step=step)
####
coef <- coef(fitlr)
theta.mat[i1-1,i2-1] <- fitlr$theta
se.theta.mat[i1-1,i2-1] <- fitlr$var.theta^.5
cor.mat[i1-1,i2-1] <- coef[1,5]
se.cor.mat[i1-1,i2-1] <- coef[1,6]
score.mat[i1-1,i2-1] <- fitlr$score
if (data.return==0) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr)
if (data.return==1) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr,data=datalr)
if (i2==2) names1 <- c(names1, paste(cut1[i1-1],"-",cut1[i1]))
if (i1==2) names2 <- c(names2, paste(cut2[i2-1],"-",cut2[i2]))
theta <- c(theta,fitlr$theta)

if ((silent<=-1) & (iid==1)) print(head(fitlr$theta.iid)); 
if ((silent<=-1) & (iid==1)) {
print(idi) ; print(datalr$tsid)
print(dim(fitlr$theta.iid))
print(head(fitlr$theta.iid))
print(dim(theta.iid))
print(length( idi %in% unique(datalr$tsid)))
}
if (iid==1) theta.iid[idi %in% unique(datalr$tsid),k] <-c(fitlr$theta.iid) 
###if (iid==1) theta.iid[rownames(fitlr$theta.iid),k] <-  fitlr$theta.iid 
}

var.thetal <- NULL
if (iid==1)  var.thetal <- t(theta.iid) %*% theta.iid

colnames(score.mat) <- colnames(cor.mat) <-  colnames(se.cor.mat)  <- colnames(se.theta.mat) <- colnames(theta.mat) <- names1; 
rownames(score.mat) <- rownames(cor.mat) <-  rownames(se.cor.mat) <-  rownames(se.theta.mat) <- rownames(theta.mat) <- names2; 

ud <- list(model.fits=ud,theta=theta.mat,var.theta=se.theta.mat^2,
	   se.theta=se.theta.mat,thetal=theta,thetal.iid=theta.iid,var.thetal=var.thetal,model=model,
	   cor=cor.mat,se.cor=se.cor.mat,score=score.mat); 
class(ud)<-"pc.twostage" 
attr(ud,"var.link")<-var.link; 
attr(ud, "Type") <- model
return(ud);
} ## }}}

##' @export
piecewise.data <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,silent=1)
{ ## {{{
ud <- list()
if (missing(cut2)) cut2 <- cut1; 
nc1 <- length(cut1); nc2 <- length(cut2)
dataud <- c()

k <- 0; 
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n"); 
 datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
			status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent) 
if (silent<=-1) print(summary(datalr)); 
if (silent<=-1) print(head(datalr)); 
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
###
datalr$strata <- paste( c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),collapse=",",sep="-")
datalr$intstrata <- 
c(paste(c(cut1[i1-1],cut1[i1]),collapse=",",sep="-"),paste( c(cut2[i2-1],cut2[i2]),collapse=",",sep="-"))

if (silent<=-1) print(head(datalr)); 
dataud <- rbind(dataud,datalr)
}

return(data.frame(dataud))
} ## }}}

##' @export
summary.pc.twostage <- function(object,var.link=NULL,...)
{ ## {{{
  if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
  
  res <- list(estimates=object$theta,se=object$se.theta,cor=object$cor,se.cor=object$se.cor,
	      model=object$model,score=object$score)
  class(res) <- "summary.pc.twostage"
  attr(res,"var.link")<-attr(object,"var.link"); 
  attr(res, "Type") <- object$model
  res
} ## }}}

##' @export
print.pc.twostage <- function(x,var.link=NULL,...)
{ ## {{{
   if (!(inherits(x,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
   print( summary(x,var.link=var.link,...))
} ## }}}

##' @export
print.summary.pc.twostage <- function(x,var.link=NULL, digits=3,...)
{ ## {{{
  
  if (is.null(var.link)) { if (attr(x,"var.link")==1) vlink <- 1 else vlink <- 0; } else vlink <- var.link
  print(vlink)

  if (x$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
  if (x$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 
 
  if (max(x$score)>0.001) { cat("Score of log-likelihood for parameter estimates (too large?)\n"); print(x$score);cat("\n\n");}

  if (vlink==1) cat("log-coefficient for dependence parameter (SE) \n")  else cat("Dependence parameter (SE) \n");
  print(coefmat(x$estimate,x$se,digits=digits,...))
  cat("\n") 

  if (x$model=="plackett") {cat("Spearman Correlation (SE) \n");cor.type <- "Spearman Correlation"; }
  if (x$model=="clayton.oakes") {cat("Kendall's tau (SE) \n"); cor.type <- "Kendall's tau";}

  print(coefmat(x$cor,x$se.cor,digits,...))
  cat("\n") 
} ## }}}

##' @export
coefmat <- function(est,stderr,digits=3,...) { ## {{{
  myest <- round(10^digits*(est))/10^digits;
  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
  mysd <- round(10^digits*(stderr))/10^digits;  
  res <- matrix(paste(format(myest)," (",format(mysd),")",sep=""),ncol=ncol(est))
  dimnames(res) <- dimnames(est)
  colnames(res) <- paste("",colnames(res))
  noquote(res)
} ## }}}

##' Fits two-stage model for describing depdendence in survival data
##' using marginals that are on cox or aalen form using the twostage funcion, but
##' call is different and easier and the data manipulation  build into the function.
##' Useful in particular for family design data. 
##'
##' If clusters contain more than two times, the algoritm uses a composite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##' @examples
##' data(prt)
##' margp<- coxph(Surv(time,status==1)~factor(country),data=prt)
##' fitco<-twostage(margp,data=prt,clusters=prt$id)
##' summary(fitco)
##' 
##' des <- model.matrix(~-1+factor(zyg),data=prt); 
##' fitco<-twostage(margp,data=prt,theta.des=des,clusters=prt$id)
##' summary(fitco)
##' 
##' dfam <- simSurvFam(1000)
##' dfam <- fast.reshape(dfam,var=c("x","time","status"))
##' 
##' desfs <- function(x,num1="num1",num2="num2")
##' { 
##' pp <- (x[num1]=="m")*(x[num2]=="f")*1   ## mother-father 
##' pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1 ## mother-child
##' cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1               ## child-child
##' c(pp,pc,cc)
##' } 
##' 
##' marg <- coxph(Surv(time,status)~factor(num),data=dfam)
##' out3 <- easy.twostage(marg,data=dfam,time="time",status="status",id="id",deshelp=0,
##'                       score.method="fisher.scoring",theta.formula=desfs,
##'                       desnames=c("parent-parent","parent-child","child-cild"))
##' summary(out3)
##' @keywords survival twostage 
##' @export
##' @param margsurv model 
##' @param data data frame
##' @param score.method Scoring method
##' @param status Status at exit time
##' @param time Exit time
##' @param entry Entry time
##' @param id name of cluster variable in data frame
##' @param Nit Number of iterations
##' @param detail Detail for more output for iterations 
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.formula design for depedence, either formula or design function
##' @param desnames names for dependence parameters
##' @param deshelp if 1 then prints out some data sets that are used, on on which the design function operates
##' @param var.link Link function for variance (exp link)
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size for newton-raphson
##' @param model plackett or clayton-oakes model
##' @param marginal.surv vector of marginal survival probabilities 
##' @param strata strata for fitting 
##' @param max.clust max clusters
##' @param se.clusters clusters for iid decomposition for roubst standard errors
easy.twostage <- function(margsurv=NULL,data=sys.parent(),score.method="nlminb",
status="status",time="time",entry=NULL,id="id", Nit=60,detail=0, silent=1,weights=NULL, control=list(),
theta=NULL,theta.formula=NULL,desnames=NULL,deshelp=0,var.link=1,iid=1,
step=0.5,model="plackett",marginal.surv=NULL,strata=NULL,max.clust=NULL,se.clusters=NULL)
{ ## {{{
### marginal trunction probabilty, to be computed from model 
pentry <- NULL

if (is.null(marginal.surv))
if (class(margsurv)[1]=="coxph")
{ ## {{{
###    ps <- survfit(margsurv)$surv
    coxformula <- margsurv$formula
    X <- model.matrix(coxformula,data=data)[,-1]; 
    baseout <- survival::basehaz(margsurv,centered=FALSE); 
    baseout <- cbind(baseout$time,baseout$hazard)
    cumh <-  Cpred(baseout,data[,time])[,2]
    RR<-exp(X %*% coef(margsurv))
    ps<-exp(-cumh*RR)
    ## }}}
  } else if (class(margsurv)[1]=="phreg")
  {  ## {{{
	  ps <- predict(margsurv)
	  pentry <- predict(margsurv,pentry)
  } ## }}} 
  else stop("marginal survival probabilities must be given as marginal.sur or margsurv \n"); 

  data <- cbind(data,ps)
  if (!is.null(pentry)) data <- cbind(data,pentry)

  ### make all pairs in the families,
  fam <- familycluster.index(data[,id])
  data.fam <- data[fam$familypairindex,]
  data.fam$subfam <- fam$subfamilyindex

  ### make dependency design using wide format for all pairs 
  data.fam.clust <- fast.reshape(data.fam,id="subfam")
  if (is.function(theta.formula)) {
     desfunction <- compiler::cmpfun(theta.formula)
    if (deshelp==1){
 	  cat("These names appear in wide version of pairs for dependence \n")
	  cat("design function must be defined in terms of these: \n")
	  cat(names(data.fam.clust)); cat("\n")
	  cat("Here is head of wide version with pairs\n")
	  print(head(data.fam.clust)); cat("\n")
    }
    des.theta  <- t( apply(data.fam.clust,1,desfunction)) 
    colnames(des.theta) <- desnames
    desnames <- desnames
     } else {
	  if (is.null(theta.formula)) theta.formula <- ~+1
          des.theta <- model.matrix(theta.formula,data=data.fam.clust)
          desnames <- colnames(des.theta); 
     }
     data.fam.clust <- cbind(data.fam.clust,des.theta)
     if (deshelp==1) {
	 cat("These names appear in wide version of pairs for dependence \n")
	     print(head(data.fam.clust))
     }

    ### back to long format keeping only needed variables
     if (is.null(pentry))
    data.fam <- fast.reshape(data.fam.clust,varying=c(id,"ps",status))
    else data.fam <- fast.reshape(data.fam.clust,varying=c(id,"ps",status,"pentry"))
    if (deshelp==1) {
	cat("Back to long format for twostage (head)\n"); 
        print(head(data.fam)); 
	cat("\n")
###	cat(paste("twostage, called with reponse",response,"\n")); 
	cat(paste("cluster=",id,",  subcluster (pairs)=subfam \n")); 
	cat(paste("design variables =")); 
	cat(desnames)
	cat("\n")
    } 

###    print(status)
###    print(names(data.fam))
###    print(data.fam[,status])
    if (is.null(pentry)) ptrunc <- NULL else ptrunc <- data.fam[,pentry]

    out <- twostage(NULL,data=data.fam,
                    clusters=data.fam$subfam,
		    theta.des=data.fam[,desnames],
                    detail=detail, score.method=score.method, Nit=Nit,step=step,
                    iid=iid,theta=theta, var.link=var.link,model=model, 
                    max.clust=max.clust,
                    marginal.survival=data.fam[,"ps"],
                    marginal.status=data.fam[,status],
		    marginal.trunc=ptrunc,
		    se.clusters=data.fam[,id])

   return(out)
} ## }}}

### library(mets)
### dfam <- simSurvFam(10000)
### dfam <- fast.reshape(dfam,var=c("x","time","status"))
###### 
### desfs <- function(x,num1="num1",num2="num2")
### { 
### pp <- (x[num1]=="m")*(x[num2]=="f")*1   ## mother-father 
### pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1 ## mother-child
### cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1               ## child-child
### c(pp,pc,cc)
### } 
######
### marg <- coxph(Surv(time,status)~factor(num),data=dfam)
### system.time(
###out3 <- easy.twostage(marg,data=dfam,time="time",status="status",id="id",deshelp=1,
###  score.method="fisher.scoring",theta.formula=desfs,
###  desnames=c("parent-parent","parent-child","child-cild"))
###)

##' @export
simSurvFam <- function(n,beta=0.0,theta=1,lam0=0.5,lam1=1,lam2=1,ctime=10,...) { ## {{{ 
###	n=10; beta=0; theta=1; lam1=1;lam2=1; ctime=10; lam0=0.5
xm <- rbinom(n,1,0.5); xf <- rbinom(n,1,0.5); 
xb1 <- rbinom(n,1,0.5); xb2 <- rbinom(n,1,0.5); 
###
zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
tm <- rexp(n)/(zf*exp(xm*beta)*lam0)
tf <- rexp(n)/(zf*exp(xf*beta)*lam0)
tb1 <- rexp(n)/((zf+zb)*exp(xb1*beta)*2*lam0)
tb2 <- rexp(n)/((zf+zb)*exp(xb2*beta)*2*lam0)
cm <- ifelse(tm<ctime,1,0); cf <- ifelse(tf<ctime,1,0); 
cb1 <- ifelse(tb1<ctime,1,0); cb2 <- ifelse(tb2<ctime,1,0); 
tm <- ifelse(tm<ctime,tm,ctime); tf <- ifelse(tf<ctime,tf,ctime)
tb1 <- ifelse(tb1<ctime,tb1,ctime); tb2 <- ifelse(tb2<ctime,tb2,ctime)
#
data.frame(xm=xm,xf=xf,xb1=xb1,xb2=xb2,timem=tm,timef=tf,timeb1=tb1,timeb2=tb2,statusm=cm,statusf=cf,
	   statusb1=cb1,statusb2=cb2,id=1:n)
} ## }}} 
