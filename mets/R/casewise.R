##' Estimates the casewise concordance based on Concordance and marginal estimate using timereg and performs test for independence
##'
##' @title Estimates the casewise concordance based on Concordance and marginal estimate using timereg and performs test for independence
##' @details Uses cluster based conservative standard errors for marginal
##' @param conc Concordance 
##' @param marg Marginal estimate
##' @param test Type of test for independence assumption. "conc" makes test on concordance scale and "case" means a test on the casewise concordance
##' @param p check that marginal probability is greater at some point than p
##' @author Thomas Scheike
##' @aliases casewise.test slope.process casewise.bin
##' @examples
##' \donttest{ ## Reduce Ex.Timings
##' data(prt);
##' 
##' prt <- prt[which(prt$id %in% sample(unique(prt$id),7500)),]
##' ### marginal cumulative incidence of prostate cancer
##' times <- seq(60,100,by=2)
##' outm <- comp.risk(Event(time,status)~+1,data=prt,cause=2,times=times)
##' 
##' cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
##' cifdz <- predict(outm,X=1,uniform=0,resample.iid=1)
##' 
##' ### concordance for MZ and DZ twins
##' cc <- bicomprisk(Event(time,status)~strata(zyg)+id(id),
##'                  data=prt,cause=c(2,2))
##' cdz <- cc$model$"DZ"
##' cmz <- cc$model$"MZ"
##' 
##' ### To compute casewise cluster argument must be passed on,
##' ###  here with a max of 100 to limit comp-time
##' outm <-comp.risk(Event(time,status)~+1,data=prt,
##'                  cause=2,times=times,max.clust=100)
##' cifmz <-predict(outm,X=1,uniform=0,resample.iid=1)
##' cc <-bicomprisk(Event(time,status)~strata(zyg)+id(id),data=prt,
##'                 cause=c(2,2),se.clusters=outm$clusters)
##' cdz <- cc$model$"DZ"
##' cmz <- cc$model$"MZ"
##' 
##' cdz <- casewise.test(cdz,cifmz,test="case") ## test based on casewise
##' cmz <- casewise.test(cmz,cifmz,test="conc") ## based on concordance
##' 
##' plot(cmz,ylim=c(0,0.7),xlim=c(60,100))
##' par(new=TRUE)
##' plot(cdz,ylim=c(0,0.7),xlim=c(60,100))
##' 
##' slope.process(cdz$casewise[,1],cdz$casewise[,2],iid=cdz$casewise.iid)
##' 
##' slope.process(cmz$casewise[,1],cmz$casewise[,2],iid=cmz$casewise.iid)
##'
##' }
##' @export
casewise.test <- function(conc,marg,test="no-test",p=0.01)
{ ## {{{
###	conc=cdz; marg=cifdz; p=0.01
###	conc=cmz; marg=cifmz
###	names(cdz)
###	cdz$casewise
  if (sum(marg$P1>p)==0) stop("No timepoints where marginal > ",p,"\n"); 
  time1 <- conc$time; time2 <- marg$time[marg$P1>0.01]
  mintime <- max(time1[1],time2[1])
  maxtime <- min(max(time1),max(time2))
  timer <- seq(mintime, maxtime,length=100)
  dtimer <- timer[2]-timer[1]

  margtime <- Cpred(cbind(marg$time,c(marg$P1)),timer)[,2]
  concP1 <-   Cpred(cbind(conc$time,c(conc$P1)),timer)[,2]
  se.margtime <- Cpred(cbind(marg$time,c(marg$se.P1)),timer)
  se.concP1 <-   Cpred(cbind(conc$time,c(conc$se.P1)),timer)
  outtest <- NULL
  casewise.iid <- NULL
  casewise <- concP1/margtime

  se.casewise2 <- as.matrix(se.concP1[,2]/margtime,ncol=100)
  se.margtime <- as.matrix(se.margtime[,2],ncol=100)
  se.casewise <- NULL;
  outtest <- NULL; 

  if (is.null(conc$P1.iid) || is.null(marg$P1.iid))   {
    cat("Warning, need iid decomposition for correct se's for Concordance \n");
  } else {
   conciid <- Cpred(cbind(conc$time,conc$P1.iid),timer)[,-1]
   nconc <- colnames(conc$P1.iid)
   P1iid <- Cpred(cbind(marg$time,marg$P1.iid),timer)[,-1]
   P1iid2  <- 2*P1iid*margtime;
   se.p2 <- apply(P1iid2^2,1,sum)^.5
   conciid <- (P1iid2[,nconc]-conciid) ### iid of conc-p1^2 test
   c.iid <- casewise.iid <- conciid/margtime - P1iid[,nconc]*casewise
   se.casewise <- apply(casewise.iid^2,1,sum)^.5
   dif.casewise.iid <- conciid/margtime  ### iid of casewise test

   if (test=="case") 
   {
        diff.iidcase <- apply(dif.casewise.iid,2,sum)*dtimer;
        sd.pepem <- sum(diff.iidcase^2)^.5
        diff.pepem <- sum((casewise-margtime))*dtimer
        z.pepem <- diff.pepem/sd.pepem
        pval.pepem <- 2*pnorm(-abs(z.pepem))
        outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem)

	### test for constant casewise concordance
	diff.const <- (casewise-mean(casewise))
	iid.constant <-  (c.iid - matrix(apply(c.iid,2,mean),nrow(c.iid),ncol(c.iid),byrow=T))
	se.const <- apply(iid.constant^2,1,sum)^.5
        test.constant <- max(abs(diff.const/se.const))
	sim.maxs <- apply(abs(iid.constant/se.const),1,max)
	pval.const <- pval(sim.maxs,test.constant)
        outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem,test.constant,pval.const)

        colnames(outtest) <- c("cum dif.","sd","z","pval","constant-case","pval") 			 
        rownames(outtest) <- "pepe-mori"

    } else if (test=="conc") {
        diff.iid <- apply(conciid,2,sum)*dtimer
        sd.pepem <- sum(diff.iid^2)^.5
        diff.pepem <- sum((concP1-margtime^2))*dtimer
        z.pepem <- diff.pepem/sd.pepem
        pval.pepem <- 2*pnorm(-abs(z.pepem))
        outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem)
        colnames(outtest) <- c("cum dif.","sd","z","pval") 			 
        rownames(outtest) <- "pepe-mori"
    }
  }

   concout <- cbind(timer,concP1,se.concP1[,2])
   colnames(concout) <- c("time","concordance","se")
   margout <- cbind(timer,margtime,se.margtime)
   colnames(margout) <- c("time","cif","se")
   casewiseout <- cbind(timer,casewise,se.casewise,se.casewise2)
   colnames(casewiseout) <- c("time","casewise","se","se2")
   difout <- cbind(timer,concP1-margtime^2,apply(conciid^2,1,sum)^.5)

  out <- list(casewise=casewiseout,marg=margout,conc=concout,casewise.iid=casewise.iid, 
	      test=outtest,mintime=mintime,maxtime=maxtime,same.cluster=TRUE,testtype=test)
  class(out) <- "casewise"
  return(out)
} ## }}}

##' @export
slope.process  <-  function(time,y,iid=NULL)
{ ## {{{ 
ctime <- scale(time)
caselm <- lm(y ~ ctime)
slope <- coef(caselm)[2]

if (!is.null(iid)) {
diff.iid <- iid 
lm.iid <- c()
for (i in 1:ncol(iid))
{
lmiid <- lm(iid[,i] ~ ctime)
lm.iid <- rbind(lm.iid,coef(lmiid))
}
se.slope <- apply(lm.iid,2,sd)^.5
} else {se.slop <- NULL}

z.slope <- slope/se.slope[2]
pval <- 2*(1-pnorm(abs(z.slope)))

out <- list(intercept=coef(caselm)[1],slope=slope,se.slope=se.slope,pval.slope=pval)
return(out)
} ## }}} 

##' .. content for description (no empty lines) ..
##'
##' @title Estimates the casewise concordance based on Concordance and marginal estimate using prodlim but no testing
##' @param conc Concordance 
##' @param marg Marginal estimate
##' @param cause.marg specififes which cause that should be used for marginal cif based on prodlim
##' @author Thomas Scheike
##' @examples
##' \donttest{ ## Reduce Ex.Timings
##' library(prodlim)
##' data(prt);
##' 
##' ### marginal cumulative incidence of prostate cancer##' 
##' outm <- prodlim(Hist(time,status)~+1,data=prt)
##' 
##' times <- 60:100
##' cifmz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="MZ")) ## cause is 2 (second cause) 
##' cifdz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="DZ"))
##' 
##' ### concordance for MZ and DZ twins
##' cc <- bicomprisk(Event(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),prodlim=TRUE)
##' cdz <- cc$model$"DZ"
##' cmz <- cc$model$"MZ"
##' 
##' cdz <- casewise(cdz,outm,cause.marg=2) 
##' cmz <- casewise(cmz,outm,cause.marg=2)
##' 
##' plot(cmz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE,col=c(3,2,1))
##' par(new=TRUE)
##' plot(cdz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE)
##' summary(cdz)
##' summary(cmz)
##' }
##' @export
casewise <- function(conc,marg,cause.marg)
{ ## {{{
  if (missing(cause.marg)) stop("Please specify cause of marginal (as given in Event object)")
  if ((!class(conc)=="prodlim")  || (!class(marg)=="prodlim")) stop("Assumes that both models are based on prodlim function \n"); 
  time1 <- conc$time
  time2 <- marg$time

  cause.prodlim <- match(as.character(cause.marg),levels(prodlim::getEvent(marg$model.response)))
  if (is.na(cause.prodlim)) stop("Cause did not match marginal model")
  
  mintime <- max(time1[1],time2[1])
  maxtime <- min(max(time1),max(time2))
  timer <- seq(mintime, maxtime,length=100)
  dtimer <- timer[2]-timer[1]
 
  out <- conc
  out$time <- timer
  if (class(marg)=="comp.risk") margtime <- Cpred(cbind(marg$time,c(marg$P1)),timer)[,2] else if (class(marg)=="prodlim") {
	  cuminc <- data.frame(marg$cuminc)[,cause.prodlim]; 
	  se.cuminc <- data.frame(marg$se.cuminc)[,cause.prodlim]; 
	  margtime <- Cpred(cbind(marg$time,c(cuminc)),timer)[,2]; 
	  se.margtime <- Cpred(cbind(marg$time,c(se.cuminc)),timer)[,2]; 
  } else stop("marginal cumulative incidence comp.risk or prodlim output\n"); 

  if (class(conc)=="comprisk") concP1 <-  Cpred(cbind(conc$time,c(conc$P1)),timer)[,2]
  else if (class(conc)=="prodlim")  {
	  conc.cuminc <- data.frame(conc$cuminc)[,1]
	  conc.se.cuminc <- data.frame(conc$se.cuminc)[,1]
          se.P1 <-  Cpred(cbind(conc$time,conc.se.cuminc),timer)[,2]
	  concP1 <- Cpred(cbind(conc$time,conc.cuminc),timer)[,2]
  }
  concordance<- cbind(timer,concP1,se.P1)
  colnames(concordance) <- c("time","cif","se.cif")

  P1 <- concP1/margtime
  se.P1 <- se.P1/margtime
  med <- (margtime>0)  & (concP1 > 0) 
  out$P1 <- P1[med]
  out$se.P1 <- se.P1[med]
  out$timer <- timer[med]


  margout <- cbind(timer,margtime,se.margtime)
  colnames(margout) <- c("time","cif","se.cif")

  probout <- cbind(out$timer,out$P1,out$se.P1)
  colnames(probout) <- c("time","casewise conc","se casewise")

  out <- list(casewise=probout,marg=margout,concordance=concordance,test=NULL)
  class(out) <- "casewise"
  return(out)
} ## }}}

##' @export
casewise.bin <- function(nc,nd)
{ ## {{{
ud <- 	glm(cbind(nc,round(0.5*nd+nc))~ +1,family=binomial())
udci <- confint(ud)
pud <- predict(ud,se.fit=TRUE,type="response")
return(list(p.casewise=pud$fit,ci.casewise=exp(udci)))
} ## }}}


##' @export
plot.casewise <- function(x,ci=NULL,lty=NULL,ylim=NULL,col=NULL,xlab="time",ylab="concordance",legend=FALSE,...)
{ ## {{{
  if (is.null(col)) col <- 1:3
  if (is.null(lty)) lty <- 1:3
  if (is.null(ylim)) ylim=range(c(x$casewise[,2],x$marg[,2]))

  plot(x$casewise[,1],x$casewise[,2],type="s",ylim=ylim,lty=lty[1],col=col[1],xlab=xlab,ylab=ylab,...)
  if (!is.null(ci)) {
     ul <- x$casewise[,2]+qnorm(1-(1-ci)/2)* x$casewise[,3]
     nl <- x$casewise[,2]-qnorm(1-(1-ci)/2)* x$casewise[,3]
     lines(x$casewise[,1],ul,type="s",ylim=ylim,lty=lty[3],col=col[3])
     lines(x$casewise[,1],nl,type="s",ylim=ylim,lty=lty[3],col=col[3])
  }
  lines(x$marg[,1],x$marg[,2],lty=lty[2],col=col[2],type="s")

  if (legend==TRUE) legend("topleft",lty=lty[1:2],col=col[1:2],c("Casewise concordance","Marginal estimate"))
} ## }}}

##' @export
summary.casewise <- function(object,marg=FALSE,...)
{ ## {{{
   cat("Casewise concordance and standard errors \n"); 
   print(signif(cbind(object$casewise),3))
   cat("\n"); 

   if (marg==TRUE) {
      cat("Marginal cumulative incidence and standard errors \n"); 
      print(signif(cbind(object$marg),3))
   }

   if (!is.null(object$test)) {
	   printcasewisetest(object)
   }

} ## }}}

##' .. content for description (no empty lines) ..
##'
##' @title prints Concordance test 
##' @param x output from casewise.test 
##' @param digits number of digits
##' @param \dots Additional arguments to lower level functions
##' @author Thomas Scheike
##' @export
printcasewisetest <- function(x,digits=3,...)
{ ## {{{
  cat("\n")
  cat("Pepe-Mori type test for H_0: conc_1(t)= conc_2(t)\n")
  if (x$same.cluster==TRUE) cat("Assuming same clusters for the two functions\n") else 
  cat("Assuming independence for estimators\n");
  cat(paste("Time.range =",signif(x$mintime,3),"--",signif(x$maxtime,3),"\n\n"));
  prmatrix(signif(x$test,digits))
  invisible(x)
}	## }}}

##' .. content for description (no empty lines) ..
##'
##' @title Concordance test Compares two concordance estimates
##' @param conc1 Concordance estimate of group 1
##' @param conc2 Concordance estimate of group 2
##' @param same.cluster if FALSE then groups are independent, otherwise estimates are based on same data. 
##' @author Thomas Scheike
##' @export
test.conc <- function(conc1,conc2,same.cluster=FALSE)
{ ## {{{
  time  <- time1 <- conc1$time
  time2 <- conc2$time
  mintime <- max(time1[1],time2[1])
  maxtime <- min(max(time1),max(time2))
  timer <- seq(mintime, maxtime,length=100)
  dtimer <- timer[2]-timer[1]

  conc2timer <- Cpred(cbind(conc2$time,c(conc2$P1)),timer)[,2]
  conc1timer <- Cpred(cbind(conc1$time,c(conc1$P1)),timer)[,2]
  outtest <- NULL

  if (is.null(conc1$P1.iid) || is.null(conc2$P1.iid)) 
	  stop("Must give iid represenation for both estimators\n");  

    if (!is.null(conc1$P1.iid) && !is.null(conc2$P1.iid))  {
    if ( ((ncol(conc1$P1.iid[,])-ncol(conc2$P1.iid[,]))!=0) && same.cluster==TRUE)
       cat("Warning, not same number of iid residuals for concordance and marginal estimate\n"); 
    }

    if (!is.null(conc1$P1.iid)) if (!is.null(conc2$P1.iid)) {
    ### iid version af integraler
      conc2P1.iid  <- Cpred(cbind(conc2$time,conc2$P1.iid[,]),timer)[,-1]
      conc1P1.iid  <- Cpred(cbind(conc1$time,conc1$P1.iid[,]),timer)[,-1]
    if ( (ncol(conc1$P1.iid)==ncol(conc2$P1.iid)) && same.cluster==TRUE) {
	diff.iid <- conc1P1.iid-conc2P1.iid
	sdiff.iid <- apply(diff.iid,2,sum)*dtimer
	sd.pepem <- sum(sdiff.iid^2)^.5
    } else {
	diff2.iid <- conc2P1.iid
	sdiff2.iid <- apply(diff2.iid,2,sum)*dtimer
	var2.pepem <- sum(sdiff2.iid^2)
	diff1.iid <- conc1P1.iid
	sdiff1.iid <- apply(diff1.iid,2,sum)*dtimer
	var1.pepem <- sum(sdiff1.iid^2)
	sd.pepem <- (var1.pepem+var2.pepem)^.5
    }
      diff.pepem <- sum(conc2timer-conc1timer)*dtimer
###	print(cbind(conc2timer,conc1timer))
      z.pepem <- diff.pepem/sd.pepem
      pval.pepem <- 2*pnorm(-abs(z.pepem))

      outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem)
      colnames(outtest) <- c("cum dif.","sd","z","pval") 			 
      rownames(outtest) <- "pepe-mori"
###  print(outtest,4)
###  prmatrix(outtest,3)			
    } 

  outtest <- list(test=outtest,mintime=mintime,maxtime=maxtime,same.cluster=same.cluster)
###attr(out,"class") <- rev(attr(out,"class")) 
  class(outtest) <- "testconc"
  return(outtest)
} ## }}}

##' convert to timereg object
##'
##' @title Convert to timereg object
##' @param obj no use
##' @author Thomas Scheike
##' @export
back2timereg <- function(obj)
{ ## {{{
  out <- obj
  attr(out,"class") <- rev(attr(out,"class")) 
  return(out)
} ## }}}

