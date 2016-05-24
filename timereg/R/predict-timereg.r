slaaop<-function(z,time,cum)
{x<-time; y<-cum; index<-sum(x<=z);if (index==0) index<-1;
retur<-y[index]; return(retur); }

pred.cum<-function(x,time,cum) {ud<-sapply(x,slaaop,time,cum); 
 return(ud)}

pred.des<-function(formula,data=sys.parent())
{ ## {{{
  call <- match.call();
  m <- match.call(expand.dots=FALSE);
  special <- c("const")
  Terms <- if(missing(data)) terms(formula, special)
           else          terms(formula, special,data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, data)
  mt <- attr(m, "terms")

  XZ<-model.matrix(Terms,m)[,drop = FALSE]
  cols<-attributes(XZ)$assign
  l.cols<-length(cols)
  semicov <- attr(Terms, "specials")$const
  ZtermsXZ<-semicov-1

  if (length(semicov)) {renaalen<-FALSE; Zterms<-c();
  for (i in ZtermsXZ) Zterms<-c(Zterms,(1:l.cols)[cols==i]);
   } else {renaalen<-TRUE;}

 if (length(semicov)) {
  X<-as.matrix(XZ[,-Zterms]);
  #covnamesX <- dimnames(XZ)[[2]][-Zterms]; dimnames(X)[[2]]<-covnamesX;
  Z<-as.matrix(XZ[,Zterms]);
  #covnamesZ <- dimnames(XZ)[[2]][Zterms];dimnames(Z)[[2]]<-covnamesZ;          
 }
  else {X<-as.matrix(XZ); #covnamesX <- dimnames(XZ)[[2]];
        Z<-FALSE; #dimnames(X)[[2]]<-covnamesX; 
  }
  X <- data.matrix(X); 
  return(list(covarX=X,covarZ=Z))
} ## }}}

aalen.des2 <-  function(formula,data=sys.parent(),model=NULL,...){ ## {{{
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$model <- NULL
  Terms <- if(missing(data)) terms(formula )
  else              terms(formula, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
 intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")

  if (model=="cox.aalen") modela <- "cox.aalen" else modela <- "aalen"
  des<-read.design(m,Terms,model=modela)
  return(des) 
} ## }}}

predict.cox.aalen <- function(object,...) predict.timereg(object,...)
predict.aalen <- function(object,...) predict.timereg(object,...)
predict.comprisk <- function(object,...) predict.timereg(object,...)

predict.timereg <-function(object,newdata=NULL,X=NULL,times=NULL,
                           Z=NULL,n.sim=500, uniform=TRUE,
                           se=TRUE,alpha=0.05,resample.iid=0,...)
{
    ## {{{
###    if (object$conv$convd>=1) stop("Model did not converge.")
 ### {{{ reading designs  and models 
  if (!(inherits(object,'comprisk') || inherits(object,'aalen')
        || inherits(object,'cox.aalen')))
        stop ("Must be output from comp.risk function")

  if(inherits(object,'aalen')) { modelType <- 'aalen';
  } else if(inherits(object,'comprisk')) { modelType <- object$model;
  } else if(inherits(object,'cox.aalen')) { 
	  if (is.null(object$prop.odds)) modelType <- 'cox.aalen' else  modelType <- 'prop.odds'; 
	  
  }
  type <- "na"
  if (modelType=="prop.odds")  type <- attr(object,'type')

  n <- length(object$B.iid) ## Number of clusters (or number of individuals
                            ## if no cluster structure is specified)

  if (se==FALSE) uniform <- FALSE
  if (is.null(object$B.iid)==TRUE & (se==TRUE | uniform==TRUE)) {
    stop("resample processes necessary for these computations, set resample.iid=1");
  }
  if (is.null(object$gamma)==TRUE) { semi<-FALSE } else { semi<-TRUE }
  ## }}} 

  ## {{{ extracts design based on the different specifications
  ## cox.aalen uses prop(...), while aalen and comp.risk use const(...)
  ## this accounts for the different number of characters
  if(inherits(object,'cox.aalen')){ indexOfFirstChar <- 6; } else { indexOfFirstChar <- 7; }
  
  ### whether or not iid time coarsening is used (only for cox-aalen)
  ### or whether time is changed due to times argument, then also changing times for iid and cum
  iidtimechange <- 0; iidtime <- 0
  if (!is.null(newdata)) {  ## {{{ newdata given 
    ##  The time-constant effects first
    formulao <- attr(object,"Formula")
    des <- aalen.des2(formula(delete.response(terms(formulao))),data=newdata,model=modelType)
    time.vars <- des$X 

    if (semi==TRUE) {
      constant.covs <- des$Z; const <- c(object$gamma)
      names(const) <-substr(dimnames(object$gamma)[[1]],indexOfFirstChar, nchar(dimnames(object$gamma)[[1]])-1)
    } else constant.covs <- NULL

    ## Then extract the time-varying effects
    ###    time.coef <- data.frame(object$cum)
    time.coef <- as.matrix(object$cum)
    if (!is.null(times)) {time.coef<-Cpred(time.coef,times); iidtimechange <- 1; iidtime <- object$cum[,1];} 
    ### SE based on iid decomposition so uses time-resolution for cox.aalen model 
    if (modelType=="cox.aalen" && (!is.null(object$time.sim.resolution)) && (se==TRUE)) 
    { iidtime <- object$time.sim.resolution; iidtimechange <- 1} 

    nobs <- nrow(newdata) 
    ## }}} 
  } else if ((is.null(Z)==FALSE) || (is.null(X)==FALSE)){ ## {{{ X, Z specified  

    if (semi) zcol <- length(c(object$gamma)) else zcol <- NULL
    if (!is.null(Z)) { prow <- nrow(Z); } 
    if (!is.null(Z)) Z <- matrix(Z,ncol=zcol)
    if (semi & is.null(Z)) Z <- matrix(0,nrow=nrow(X),ncol=zcol); 

    xcol<-ncol(object$cum)-1
    if (!is.null(X)) X <- matrix(X,ncol=xcol)
    else {
	    X <-   matrix(0,xcol,1)
	    X[,1] <- 1
    }
    if (semi & is.null(Z)) Z <- matrix(0,nrow=nrow(X),ncol=zcol); 
    time.vars <- X
    if (semi) constant.covs <- Z else constant.covs <- NULL

    nobs<-nrow(X);

    ## Then extract the time-varying effects
    time.coef <- as.matrix(object$cum)

    if (!is.null(times)) {time.coef<-Cpred(time.coef,times); iidtimechange <- 1; iidtime <- object$cum[,1];} 
    ### SE based on iid decomposition so uses time-resolution for cox.aalen model 
    if (modelType=="cox.aalen" && (!is.null(object$time.sim.resolution)) && (se==TRUE)) 
    { iidtime <- object$time.sim.resolution; iidtimechange <- 1} 

    ## }}} 
  } else { ## {{{ 
    stop("Must specify either newdata or X, Z\n");
  } ## }}} 

  ## }}}

  ## {{{ predictions for competing risks and survival data

  cumhaz<-as.matrix(time.vars) %*% t(matrix(time.coef[,-1],ncol=(ncol(time.coef)-1)))
  times <- time<-time.coef[,1]; 
  if (semi==TRUE) pg <- nrow(object$gamma); 
  nt<-length(time);

  ### set up articial time.pow for aalen  and cox.aalen to make unified code for 
  ###  comp.risk and survival
  if(inherits(object,'aalen') & semi==TRUE) timepow <- rep(1,pg)
  if(inherits(object,'cox.aalen')) timepow <- rep(0,pg)
  if(inherits(object,'comprisk')) timepow <- attr(object,"time.pow")

  if (semi==TRUE)
  constant.part <- constant.covs %*% ((matrix(rep(c(time),pg),pg,nt,byrow=TRUE)^timepow)*c(object$gamma))

  if (inherits(object,'comprisk')) { ## {{{ competing models
    if (modelType == "additive") {
      if (semi==FALSE){
        P1=1-exp(-cumhaz);
      } else {
        P1=1-exp(-cumhaz-constant.part )
      }
      RR<-1; 
    } else if (modelType == 'rcif') { # P1=exp(x^T b(t) + z^t t^p gamma) 
        if (semi==FALSE){
           P1=exp(cumhaz);
         } else {
         P1<-exp(cumhaz+constant.part);
       }
       RR<-1;
    } else if (modelType == 'rcif2') { # P1=x^T b(t) exp( z^t t^p gamma) 
        if (semi==FALSE){
         P1=cumhaz;
         RR<-1;
         } else {
         P1<-cumhaz*exp(constant.part);
         RR<-exp(constant.part);
       }
    } else if (modelType == 'prop') {# model proportional , Fine Gray extension
        if (semi==FALSE){
        RR<-exp(cumhaz);
      } else {
        RR<-exp(cumhaz+constant.part);
    }
    P1<-1-exp(-RR);
    } else if (modelType == 'fg') {# model proportional, Fine-Gray parametrization 
        if (semi==FALSE){
        RR<-cumhaz;
      } else {
        RR<-cumhaz*exp(constant.part);
    }
    P1<-1-exp(-RR);
    } else if (modelType == 'logistic') { #model logistic
      if (semi==FALSE){ RR<-exp(cumhaz); }   else { RR<-exp(cumhaz+constant.part); }
      P1<-RR/(1+RR);
    } else if (modelType == 'logistic2') { #model logistic, baseline-par
      if (semi==FALSE){ RR<-1; }   else { RR<-exp(constant.part); }
      P1<-RR*cumhaz/(1+RR*cumhaz);
    } ## }}}
    } 
    else if (modelType=="prop.odds")
    {
       RR <- exp(constant.part)
       HRR <- cumhaz* RR 
       P1 <- HRR/(1+HRR) 
       S0 <- 1/(1+HRR) 
    }
    else {  # aalen or cox.aalen survival model  ## {{{
       if (modelType == "aalen") {    #Aalen model
          if (semi==FALSE){ S0=exp(-cumhaz); } else { S0=exp(-cumhaz-constant.part) }
       RR<-NULL; 
       } else if(modelType == 'cox.aalen'){  #Cox-Aalen model
       if(semi == FALSE){ RR <- NULL; S0 <- exp(-cumhaz); } else {
         RR <- exp(constant.part);
         S0 <- exp(-cumhaz * RR);
       }
     }
     else stop("model class not supported by predict\n")
    } ## }}}

    ## }}}

  se.P1 <- NULL
  se.S0 <- NULL
  P1.iid <- NULL
  S0.iid <- NULL
  uband <- NULL
  ## i.i.d decomposition for computation of standard errors  ## {{{
  if (se==1) {
    pg<-length(object$gamma); 
    delta<-c();
    for (i in 1:n) {
       if (iidtimechange==1) 
       tmptiid<- t(Cpred(cbind(iidtime,object$B.iid[[i]]),times)[,-1,drop=FALSE])
       else tmptiid <- t(object$B.iid[[i]])
       tmp<- as.matrix(time.vars) %*% tmptiid

       if (semi==TRUE) {
             gammai <- matrix(object$gamma.iid[i,],pg,1); 
             tmp.const<-constant.covs %*% ((matrix(rep(c(time),pg),pg,nt,byrow=TRUE)^timepow)*c(gammai))
       } 

      if (i==0) { ## {{{ test print stuff
        print(tmp.const);
        if (modelType=="additive" || modelType == 'aalen'){ 
          print(tmp.const %*% matrix(time,1,nt))
        } else if (modelType=="prop"){
          print(tmp.const %*% matrix(1,1,nt));
        } else if (modelType=="cox.aalen") {
          tmp <- RR * tmp + RR * cumhaz * matrix(tmp.const,nobs,nt);
        } else if (modelType=="prop.odds") {

	}
      } ## }}} 

      if (semi==TRUE){
        if(modelType=="additive" || modelType == "aalen") { tmp<-tmp+ tmp.const } 
	else if (modelType=="prop" || modelType=="rcif") { tmp<-RR*tmp+RR*tmp.const; } 
	else if (modelType=="logistic" || modelType=="rcif2") { tmp<-RR*tmp+RR*cumhaz*tmp.const; } 
	else if (modelType=="logistic2") { tmp<-RR*tmp+RR*cumhaz*tmp.const; } 
	else if (modelType=="cox.aalen") { tmp <- RR * tmp + RR * cumhaz * tmp.const }
	else if (modelType=="prop.odds") { tmp <- RR * tmp + RR * cumhaz * tmp.const; }
      } else {
	if (modelType=="prop") { tmp<-RR*tmp; } 
      }

      delta<-cbind(delta,c(tmp)); 
    }
    se<-apply(delta^2,1,sum)^.5

    if(modelType == 'additive' || modelType == 'prop' || modelType=="fg"){ 
       se.P1<-matrix(se,nobs,nt)*(1-P1); 
       if (resample.iid==1)  P1.iid <- array(delta*c(1-P1),c(nobs,nt,n));   
    } 
    else if(modelType == 'rcif' ){ 
       se.P1<-matrix(se,nobs,nt)*P1 
       if (resample.iid==1) P1.iid <- array(delta*P1,c(nobs,nt,n));   
    } 
    else if(modelType == 'rcif2'){ 
       se.P1<-matrix(se,nobs,nt) 
       if (resample.iid==1) P1.iid <- array(delta,c(nobs,nt,n));   
    }
    else if (modelType == 'logistic'){ 
       se.P1<-matrix(se,nobs,nt)*P1/(1+RR) 
       if (resample.iid==1) P1.iid <- array(delta*c(P1/(1+RR),c(nobs,nt,n)));   
    } 
    else if (modelType == 'logistic2'){ 
       se.P1<-matrix(se,nobs,nt)*1/(1+cumhaz*RR)^2 
       if (resample.iid==1) P1.iid <- array(delta*c(1/(1+cumhaz*RR),c(nobs,nt,n)));   
    } 
    else if (modelType == 'aalen' || modelType == 'cox.aalen'){ 
       se.S0<-matrix(se,nobs,nt)*S0 
       if (resample.iid==1) S0.iid <- array(delta*c(S0),c(nobs,nt,n));   
    }
    else if (modelType == 'prop.odds'){ 
       if (attr(object,'type')=="comprisk") {
          se.P1 <-matrix(se,nobs,nt)*S0^2
          if (resample.iid==1) P1.iid <- array(delta*c(S0^2),c(nobs,nt,n));   
       }
       if (attr(object,'type')=="survival") {
          se.S0<-matrix(se,nobs,nt)*S0^2 
          if (resample.iid==1) S0.iid <- array(delta*c(S0^2),c(nobs,nt,n));   
       }
    }
    }
    ## }}}

    ### uniform confidence bands, based on resampling  ## {{{
    if (uniform==1) {
      mpt <- .C('confBandBasePredict',
                delta = as.double(delta), nObs = as.integer(nobs), nt = as.integer(nt),
                n = as.integer(n), se = as.double(se), mpt = double(n.sim*nobs),
                nSims = as.integer(n.sim), PACKAGE="timereg")$mpt;
  
      mpt <- matrix(mpt,n.sim,nobs,byrow = TRUE);
      uband <- apply(mpt,2,percen,per=1-alpha);
    } else uband<-NULL; 
   ## }}}

  if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic"
     || modelType=='rcif2' || modelType=='rcif' || modelType=='fg' || modelType=='logistic2'){
    P1<-matrix(P1,nrow=nobs);
  } else if (modelType == 'aalen' || modelType == 'cox.aalen'){
    S0<-matrix(S0,nrow=nobs);
  }
  else if (modelType == 'prop.odds'){
   if (attr(object,'type')=="comprisk") P1<-matrix(P1,nrow=nobs);
   if (attr(object,'type')=="survival") S0<-matrix(S0,nrow=nobs);
  }

  out<-list(time=time,unif.band=uband,model=modelType,alpha=alpha,
            newdata=list(X = time.vars, Z = constant.covs),RR=RR,
            call=sys.calls()[[1]], initial.call = attr(object,'Call'));

  if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic"
     || modelType=='rcif2' || modelType=='rcif' || modelType=='fg' || modelType=='logistic2' ||
     ((modelType=='prop.odds') && type=="comprisk")){
    if (nrow(P1)==1)  { P1 <- c(P1); se.P1 <- c(se.P1); }
    out$P1 <- P1;
    out$se.P1 <- se.P1;    
    out$clusters <- attr(object,"clusters"); 
    if (resample.iid==1) {out$P1.iid <- P1.iid[1,,]; colnames(out$P1.iid)<-paste(unique(out$clusters));}
  } else if (modelType == 'aalen' || modelType == 'cox.aalen' || 
     ((modelType=='prop.odds') && type=="survival")){
    out$S0 <- S0;
    out$se.S0 <- se.S0;    
    if (resample.iid==1) {out$S0.iid <- S0.iid[1,,]; colnames(out$S0.iid)<-paste(unique(out$clusters));}
  }
   # e.g. for an compound risk model, className = predictComprisk
  className <- switch(class(object),aalen='predictAalen',cox.aalen='predictCoxAalen',comprisk='predictComprisk')


  subclass <- switch(type,comprisk="comprisk",survival="survival",na="na")
  class(out) <- "predict.timereg"
  attr(out,'className') <- className
  attr(out,'subclass') <- subclass

  return(out)
} ## }}}

pava <- function(x, w=rep(1,length(x)))  # R interface to the compiled code
{ ## {{{
  n = length(x)
  if (n != length(w)) return (0)    # error
  result  = .C("pava",
        y = as.double(x),
        as.double(w),
        as.integer(n) )
  result[["y"]]
} ## }}}

plot.predict.timereg<-function(x,uniform=1,new=1,se=1,col=1,lty=1,lwd=2,multiple=0,specific.comps=0,ylim=c(0,1),
xlab="Time",ylab="Probability",transparency=FALSE,monotone=TRUE,...)
{ ## {{{
  object <- x; rm(x);
  modelType <- object$model
  time<-object$time;
  uband<-object$unif.band;
  nobs<-nrow(object$newdata$X);
  RR<-object$RR;
  alpha <- object$alpha;
  ### Here we use mainLine as the central line (between confidence
  ### intervals or bands), so that we don't have to distinguish
  ### between the case when we want to plot a predicted survival function
  ### and the case when we want to plot a predicted risk funcion
  
  subtype <- attr(object,'subclass')

###  if  ((modelType=='prop.odds')) {
###	  subtype <- attr(object,'type'); 
###  } else subtype <- ""
###  print(modelType)
###  print(subtype) 

  if (modelType == 'aalen' || modelType == 'cox.aalen' ||
     ((modelType=='prop.odds') && subtype=='survival')){
    type<-"surv"
    mainLine <- as.matrix(object$S0);
    if (monotone==TRUE) { mainLine<--t(apply(as.matrix(-mainLine),1,pava)); 
    mainLine[mainLine<0]<-0; 
    mainLine[mainLine>1]<-1; 
    }
    if (is.null(object$se.S0))  mainLine.se <- NULL else mainLine.se <- as.matrix(object$se.S0);    
  } else if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic"
     || modelType=='rcif2' || modelType=='rcif' || modelType=='fg' || modelType=='logistic2' ||
     ((modelType=='prop.odds') && subtype=='comprisk')){
    type<-"cif"
    mainLine <- as.matrix(object$P1);
    if (monotone==TRUE) { mainLine<-t(apply(as.matrix(mainLine),1,pava)); 
                           mainLine[mainLine<0]<-0; 
                           mainLine[mainLine>1]<-1; 
    }
    if (is.null(object$se.P1))  mainLine.se <- NULL else {
    mainLine.se <-matrix(object$se.P1,nrow=nrow(mainLine));    
    }
  }
###  print(head(mainLine))
###  print(dim(mainLine))
###  print(object$se.P1)
###  print(head(mainLine.se))
###  print(dim(mainLine.se))

  if (length(col)!=nobs){ col<-rep(col[1],nobs); }
  if (length(lty)!=nobs){ lty<-rep(lty[1],nobs); }
  if (length(lwd)!=nobs){ lwd<-rep(lwd[1],nobs); }
  if (length(uniform)!=nobs){ uniform<-rep(uniform[1],nobs); }
  if (length(se)!=nobs){ se <-rep(se[1],nobs); }
  if (sum(specific.comps)==0){
    comps<-1:nobs
  } else {
    comps<-specific.comps
  }

  for (i in comps) {
    if (new==1 & (multiple!=1 | i==comps[1])) {
      plot(time,mainLine[i,],type="s",xlab=xlab,ylab=ylab,col=col[i],
		      lty=lty[i],lwd=lwd[i],ylim=ylim,...)
    } else {
      lines(time,mainLine[i,],type="s",col=col[i],lty=lty[i],lwd=lwd[i])
    }

    if (se[1]>=1 & is.null(mainLine.se)==FALSE ) {
      lower<-mainLine[i,]-qnorm(1-alpha/2)*mainLine.se[i,]
      upper<-mainLine[i,]+qnorm(1-alpha/2)*mainLine.se[i,]
       if (monotone==TRUE) { 
       if (type=="cif") { lower<- pava(lower); upper<- pava(upper); }
       if (type=="surv") { lower<- -pava(-lower); upper<- -pava(-upper); }
        lower[lower<0]<-0; lower[lower>1]<-1; 
        upper[upper<0]<-0; upper[upper>1]<-1; 
       }

      lines(time,lower,type="s",col=col[i],lty=se[i],lwd=lwd[i]/2);
      lines(time,upper,type="s",col=col[i],lty=se[i],lwd=lwd[i]/2);
    }

    if (uniform[1]>=1 & is.null(uband)==FALSE ) {
      #if (level!=0.05) c.alpha<-percen(object$sim.test[,i],1-level)
      #else c.alpha<-object$conf.band.cumz[i];
      c.alpha=uband[i]; 
      upper<-mainLine[i,]-uband[i]*mainLine.se[i,];
      lower<-mainLine[i,]+uband[i]*mainLine.se[i,];
       if (monotone==TRUE) { 
          if (type=="cif") { lower<- pava(lower); upper<- pava(upper); }
          if (type=="surv") { lower<- -pava(-lower); upper<- -pava(-upper); }
          lower[lower<0]<-0; lower[lower>1]<-1; 
          upper[upper<0]<-0; upper[upper>1]<-1; 
       }
      if (transparency==0 || transparency==2) {
      lines(time,upper,type="s",col=col[i],lty=uniform[i],lwd=lwd[i]/2);
      lines(time,lower,type="s",col=col[i],lty=uniform[i],lwd=lwd[i]/2);
      }

    ## Prediction polygons bandds ## {{{
    if (transparency>=1) {
     col.alpha<-0.2
     col.ci<-"darkblue"
     col.ci<-col[i]; 
     lty.ci<-2
      if (col.alpha==0) col.trans <- col.ci
      else
      col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))

      #print(t); print(ci)
      n<-length(time)
      tt<-seq(time[1],time[n],length=n*10); 
      ud<-Cpred(cbind(time,upper,lower),tt)[,2:3]
      tt <- c(tt, rev(tt))
      yy <- c(upper, rev(lower))
#      tt <- c(time, rev(time))
#      yy <- c(upper, rev(lower))
     yy <- c(ud[,1], rev(ud[,2]))
      polygon(tt,yy, col=col.trans, lty=0)      
  } ## }}}

    }
  }
} ## }}}

print.predict.timereg <- function(x,...){ ## {{{

  object <- x; rm(x);
  if(!(inherits(object,'predict.timereg') )) stop('Wrong class of object');
###	  || inherits(object,'predictCoxAalen') ||
###       inherits(object,'predictComprisk'))){
###    stop('Wrong class of object');
###  }

  if (is.null(object$newdata$Z)==TRUE) semi<-FALSE else semi<-TRUE
  
  modelType <- object$model;
  modelAnnouncement <- ' Predicted survival for'
  addTo <- switch(modelType,
                  cox.aalen = 'a Cox-Aalen hazard model',
                  aalen     = 'an Aalen hazard model',
                  prop = 'a proportional competing risks (Fine-Gray type)',
                  fg    =  'a proportional competing risks (Fine-Gray type)',
                  rcif  = 'a proportional risk competing risks',
                  rcif2 = 'a proportional risk competing risks',
                  logistic = 'a logistic competing risks',
                  logistic2 = 'a logistic competing risks',
                  additive = 'an additive competing risks')
  modelAnnouncement <- paste(modelAnnouncement,addTo,'model',sep = ' ')
  cat(modelAnnouncement,fill=TRUE)

  cat(" Nonparametric terms : "); cat(colnames(object$newdata$X)[-1]); cat("   \n");  
  if (semi == TRUE) {
    cat(" Parametric terms :  "); cat(rownames(object$newdata$Z)); 
    cat("   \n");  } 
  cat("   \n");  
  
  call <- object$call;
  cat('Call to predict:',fill=TRUE);
  print(call)
  call <- object$initial.call;
  cat('Initial call:',fill=TRUE);
  print(call)
    
} ## }}}

summary.predict.timereg <- function(object,...){ ## {{{

  if(!(inherits(object,'predict.timereg') )) stop('Wrong class of object');
###  if(!(inherits(object,'predictAalen') ||
###       inherits(object,'predictCoxAalen') ||
###       inherits(object,'predictComprisk'))){
###    stop('Wrong class of object');
###  }

  modelClass <- class(object)
  modelType <- object$model;
  time<-object$time;
  uband<-object$unif.band;
  nobs<-nrow(object$newdata$X);
  RR<-object$RR;
  alpha <- object$alpha;
  call <- object$call;
  if (modelType == 'aalen' || modelType == 'cox.aalen'){
    se <- object$se.S0;    
  } else if(modelType == 'additive' || modelType == 'prop'){
    se <- object$se.P1;    
  }
  
  modelAnnouncement <- 'Predicted survival for'
    addTo <- switch(modelType,
                  cox.aalen = 'a Cox-Aalen hazard model',
                  aalen     = 'an Aalen hazard model',
                  prop = 'a proportional competing risks (Fine-Gray type)',
                  fg    =  'a proportional competing risks (Fine-Gray type)',
                  rcif  = 'a proportional risk competing risks',
                  rcif2 = 'a proportional risk competing risks',
                  logistic = 'a logistic competing risks',
                  logistic2 = 'a logistic competing risks',
                  additive = 'an additive competing risks')
  modelAnnouncement <- paste(modelAnnouncement,addTo,'model',sep = ' ')
  cat(modelAnnouncement,fill=TRUE)
  timeStatement <- paste('At',length(time),'times:',paste(c(head(time),''),collapse = ', '),'...,',time[length(time)])
  cat(timeStatement,fill=TRUE)
  obsStatement <- paste('Given covariates for',nobs,'new observations');
  cat(obsStatement,fill=TRUE)
  if(is.null(se)){
    addTo <- " - not yet done";
  } else if(is.null(uband)){
    addTo <- " - only pointwise calculations have been done"
  } else {
    addTo <- " - pointwise CI and uniform confidence band available"
  }
  cat('Standard error calculations:',fill=TRUE);
  cat(addTo,fill=TRUE);
  cat('Call:',fill=TRUE);
  print(call)
  
} ## }}}

