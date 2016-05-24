##' Estimation of twostage model with cluster truncation in bivariate situation 
##'
##' @title Estimation of twostage model with cluster truncation in bivariate situation 
##' @param survformula Formula with survival model aalen or cox.aalen, some limitiation on model specification due to call of fast.reshape (so for example interactions and * and : do not work here, expand prior to call)
##' @param data Data frame
##' @param theta.des design for dependence parameters in two-stage model 
##' @param clusters clustering variable for twins
##' @param Nit number of iteration 
##' @param final.fitting TRUE to do final estimation with SE and ... arguments for marginal models 
##' @param ... Additional arguments to lower level functions
##' @author Thomas Scheike
##' @export
##' @examples
##' data(diabetes)
##' v <- diabetes$time*runif(nrow(diabetes))*rbinom(nrow(diabetes),1,0.5)
##' diabetes$v <- v
##' 
##' aout <- twin.clustertrunc(Surv(v,time,status)~1+treat+adult,
##'			 data=diabetes,clusters="id")
##' aout$two        ## twostage output
##' par(mfrow=c(2,2))
##' plot(aout$marg) ## marginal model output
##' 
##' out <- twin.clustertrunc(Surv(v,time,status)~1+prop(treat)+prop(adult),
##'			 data=diabetes,clusters="id")
##' out$two        ## twostage output
##' plot(out$marg) ## marginal model output
twin.clustertrunc <- function(survformula,data=sys.parent(),theta.des=NULL,clusters=NULL,
		      Nit=10,final.fitting=FALSE,...)
{ ## {{{ 
## {{{  adding names of covariates from survival model to data frame if needed
## adds names that are not in data (typically intercept from additive) or 
### expansion of factors, 
## also reducing only to needed covariates 

survnames <-  all.vars(update(survformula, .~0)) 
if (length(survnames)!=3) stop("Must give entry, exit and status")
entry <- survnames[1]
exit <- survnames[2]
status <- survnames[3]


tss <- terms(survformula)
Znames <- attr(tss,"term.labels")
ZnamesO <- Znames
### checks if cluster is given in survformula then removes 
clustervar <- grep("^cluster[(][A-z0-9._:]*",Znames,perl = TRUE)
if (length(clustervar)>=1) Znames <- Znames[-clustervar]
propvar <- grep("^prop[(][A-z0-9._:]*",Znames,perl = TRUE)
if (length(propvar)>=1) model <- "cox.aalen" else model <- "aalen"

### removes prop from names
if (model=="cox.aalen") {
   Zn <- c()
   nn <- length(Znames)
   ### droppe prop() for navne 
   for (i in 1:nn) Zn <- c(Zn,substr(Znames[i],6,nchar(Znames[i])-1))
   Znames <- Zn
} 

clust.orig <- data[,clusters]
d0 <- data[,c(entry,exit,status,clusters,Znames)]
data <- d0
data$dataid <- 1:nrow(data)

d2 <- fast.reshape(data,id=clusters)
d2 <- na.omit(d2)
### only double entry people 
data <- fast.reshape(d2,labelnum=TRUE)

des <- aalen.des(survformula,data=data,model=model)
factornamesX  <-  !(des$covnamesX %in% Znames)
colnames(des$X) <- des$covnamesX
if (sum(factornamesX)>=1) data <- cbind(data,des$X[,factornamesX,drop=FALSE])
if (model=="cox.aalen") {
	factornamesZ   <-  !(des$covnamesZ %in% Znames)
        colnames(des$Z) <- des$covnamesZ
        if (sum(factornamesZ)>=1) data <- cbind(data,des$Z[,factornamesZ,drop=FALSE])
}
namesX <- des$covnamesX
namesZ <- des$covnamesZ

pweight <- rep(1,nrow(data))
if (!is.null(clusters)) nclusters <- data[,clusters] else stop("must give clusters\n"); 
if (is.null(theta.des)) ptheta <- 0 else ptheta <- rep(0,ncol(theta.des))

###singletons might be dropped so same for theta.des
theta.des <- theta.des[data$dataid,]
## }}} 
 
  assign("pweight",pweight,envir=environment(survformula))

for (i in 1:Nit)
{ ## {{{ 
  if (model=="cox.aalen") { aout <- cox.aalen(survformula,data=data,weights=1/pweight,robust=0,n.sim=0);
                            beta <- c(aout$gamma,aout$cum[,-1])
  }  else  {
           aout <- aalen(survformula,data=data,weights=1/pweight,robust=0,n.sim=0);
           beta <- aout$cum[,-1]
  }
  if (i==1) { 
   if (model=="cox.aalen") pbeta <- 0*c(aout$gamma,aout$cum[,-1]) else pbeta <- 0*aout$cum[,-1]
  } 
  if (i>=2) tout <- two.stage(aout,data=data,clusters=nclusters,theta.des=theta.des,theta=ptheta)
  else tout <- two.stage(aout,data=data,clusters=nclusters,theta.des=theta.des)
  if (!is.null(theta.des)) theta <- c(theta.des %*% tout$theta) else theta <- tout$theta
  if (attr(tout, "var.link") == 1) theta <- exp(tout$theta) 
###
  data$thetades  <- c(theta)
  dd <- fast.reshape(data,id=clusters)
  v1 <- dd[,paste(entry,"1",sep="")]; 
  v2 <- dd[,paste(entry,"2",sep="")]
  time1 <- dd[,paste(exit,"1",sep="")]; 
  time2 <- dd[,paste(exit,"2",sep="")]
  status1 <- dd[,paste(status,"1",sep="")]; 
  status2 <- dd[,paste(status,"2",sep="")]
  if (model=="cox.aalen")
  { ## {{{  prediction using covariates of subject 1 and 2
  Z1 <- as.matrix(dd[,paste(namesZ,"1",sep="")])
  Z2 <- as.matrix(dd[,paste(namesZ,"2",sep="")])
###  print(head(Z1)); print(head(Z2))
  ptruncv1t2 <- predict.two.stage(tout,Z=Z1,Z2=Z2,times=v1,times2=time2,theta=theta)
  ptrunct1v2 <- predict.two.stage(tout,Z=Z1,Z2=Z2,times=time1,times2=v2,theta=theta)
  } else {
  X1 <- as.matrix(dd[,paste(namesX,"1",sep="")])
  X2 <- as.matrix(dd[,paste(namesX,"2",sep="")])
###  print(head(X1)); print(head(X2))
  ptruncv1t2 <- predict.two.stage(tout,X=X1,X2=X2,times=v1,times2=time2,theta=theta)
  ptrunct1v2 <- predict.two.stage(tout,X=X1,X2=X2,times=time1,times2=v2,theta=theta)
  } ## }}} 
  nn <- nrow(dd)
  ppv1t2 <- .Call("claytonoakesR",1/dd$thetades1,rep(0,nn),status2,
		  ptruncv1t2$S1,ptruncv1t2$S2)$like
  ppt1v2 <- .Call("claytonoakesR",1/dd$thetades1,status1,rep(0,nn),
		  ptrunct1v2$S1,ptrunct1v2$S2)$like
  ppt1v2[status1==0] <- ppt1v2[status1==0]/ptrunct1v2$S1[status1==0]
  ppv1t2[status2==0] <- ppv1t2[status2==0]/ptruncv1t2$S2[status2==0]
###
  dd$weight1 <- c(ppt1v2)
  dd$weight2 <- c(ppv1t2)
###  head(dd); dim(dd)
  dd2 <- fast.reshape(dd,labelnum=TRUE)
  pweight <- dd2$weight
  dtheta <- sum(abs(tout$theta-ptheta))
  dmarg <- sum(abs(beta-pbeta))
  if ((dtheta+dmarg) < 0.001) break; 
  ptheta <- tout$theta
  if (model=="aalen") pbeta <- aout$cum[,-1] else  pbeta <- c(aout$gamma,aout$cum[,-1])
###  print(c(dtheta,dmarg))
} ## }}} 

if (final.fitting==TRUE) {  ## {{{ 
  if (model=="cox.aalen") aout <- cox.aalen(survformula,data=data,weights=1/pweight,n.sim=0,...) else 
           aout <- aalen(survformula,data=data,weights=1/pweight,n.sim=0,...)
  tout <- two.stage(aout,data=data,clusters=nclusters,theta.des=theta.des)
} ## }}} 

res <- list(marg=aout,two=tout,marg.weights=pweight,dtheta=dtheta,dmarg=dmarg,model=model)

return(res)
} ## }}} 

##' @export
###twin.dobdata <- function(survformula,data=data,clusters=NULL,
###                         entry="v",exit="time",status="status")
###{ ## {{{ 
##### {{{  adding names of covariates from survival model to data frame if needed
##### adds names that are not in data (typically intercept from additive) model
##### also reducing only to needed covariates 
###
###d0 <- data[,c(entry,exit,status,clusters)]
###
###model <- "aalen"
###des <- aalen.des(survformula,data=data,model=model)
###X <- des$X
######med <- des$covnamesX %in% names(data)
###   colnames(X) <- des$covnamesX
###   d0 <- cbind(d0,X)
###
###data <- d0
###
###d2 <- fast.reshape(data,id=clusters)
###d2 <- na.omit(d2)
###data <- fast.reshape(d2,numlabel=TRUE)
###
######theta.des <- data.frame(theta.des)
######theta.des$clusters <- data[,clusters]
######
######t2 <- fast.reshape(theta.des,id=clusters)
######t2 <- na.omit(t2)
######theta.des <- fast.reshape(t2,numlabel=TRUE)
###
###return(data)
###}
###
