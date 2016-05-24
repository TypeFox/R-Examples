##' Internal function.
##' Calculates Inverse Probability of Censoring
##' Weights (IPCW) and adds them to a data.frame
##'
##' @title Inverse Probability of Censoring Weights
##' @param formula Formula specifying the censoring model 
##' @param data data frame
##' @param cluster clustering variable
##' @param same.cens For clustered data, should same censoring be assumed (bivariate probability calculated as mininum of the marginal probabilities)
##' @param obs.only Return data with uncensored observations only
##' @param weight.name Name of weight variable in the new data.frame
##' @param indi.weight Name of individual censoring weight  in the new data.frame
##' @param trunc.prob If TRUE truncation probabilities are also calculated and stored in 'weight.name2' (based on Clayton-Oakes gamma frailty model)
##' @param weight.name2 Name of truncation probabilities
##' @param cens.model Censoring model (default Aalens additive model)
##' @param pairs For paired data (e.g. twins) only the complete pairs are returned (With pairs=TRUE)
##' @param theta.formula Model for the dependence parameter in the Clayton-Oakes model (truncation only)
##' @param ... Additional arguments to censoring model 
##' @author Klaus K. Holst
##' @examples
##' \dontrun{
##' data(prt)
##' prtw <- ipw(Surv(time,status==0)~country, data=prt[sample(nrow(prt),5000),],
##'             cluster="id",weight.name="w")
##' plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
##' count <- 0
##' for (l in unique(prtw$country)) {
##'     count <- count+1
##'     prtw <- prtw[order(prtw$time),]
##'     with(subset(prtw,country==l),
##'          lines(time,w,col=count,lwd=2))
##' }
##' legend("topright",legend=unique(prtw$country),col=1:4,pch=-1,lty=1)
##' }
##' @export
ipw <- function(formula,data,cluster,
                same.cens=FALSE,obs.only=TRUE,
                weight.name="w",
                trunc.prob=FALSE,weight.name2="wt",indi.weight="pr",
                cens.model="aalen", pairs=FALSE,
                theta.formula=~1, ...) { 
## {{{
    ##iid=TRUE,

    ##cens.args <- c(list(formula,n.sim=0,robust=0,data=eval(data)),list(...))
    if (tolower(cens.model)%in%c("weibull","phreg.par","phreg.weibull")) {
	ud.cens <- phreg.par(formula,data,...)
        pr <- predict(ud.cens)
        noncens <- which(!ud.cens$status)        
    } else { ## {{{ 
        m <- match.call(expand.dots = FALSE)
        m <- m[match(c("", "formula", "data", "subset", "na.action"), 
                     names(m), nomatch = 0)]
        special <- c("strata", "factor", "NN", "cluster", "dummy")
        if (missing(data)) 
            Terms <- terms(formula)
        else Terms <- terms(formula, data = data)
        m$formula <- Terms
        m[[1]] <- as.name("model.frame")
        M <- eval(m, parent.frame())
        censtime <- model.extract(M, "response")
        if (ncol(censtime)==3) {
            status <- censtime[,3]
            otimes <- censtime[,2]
            ltimes <- censtime[,1]
        } else {
            status <- censtime[,2]
            otimes <- censtime[,1]
        }
        noncens <- !status
        if (is.null(attr(terms(formula,"prop"),"specials")$prop)) {
            ud.cens <- aalen(formula,n.sim=0,robust=0,data=data,...)
            XZ <- model.matrix(formula,data)
###         Gcx <- ud.cens$cum[prodlim::sindex(ud.cens$cum[,1],otimes),-1]
            Gcx<-Cpred(ud.cens$cum,otimes)[,-1];            
            Gcx<-exp(-apply(Gcx*XZ,1,sum))            
        } else {
            ud.cens <- cox.aalen(formula,n.sim=0,robust=0,data=data,...)
            XZ <- model.matrix(formula,data)
            Gcx<-Cpred(ud.cens$cum,otimes)[,-1];
            Gcx<-exp(-apply(Gcx*XZ,1,sum))            
        }
        ##ud.cens <- do.call(cens.model,cens.args)
        Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
        pr <- Gcx
	data[,indi.weight] <- Gcx
    } ## }}} 

     if (trunc.prob) stop("Under development\n"); 
    if (trunc.prob & ncol(censtime)==3) { ## truncation ## {{{ 
        data$truncsurv <- Surv(ltimes,otimes,noncens)
        trunc.formula <- update(formula,truncsurv~.)        
        ud.trunc <- aalen(trunc.formula,data=data,robust=0,n.sim=0,residuals=0,silent=1)
        dependX0 <- model.matrix(theta.formula,data)        
        twostage.fit <-two.stage(ud.trunc,data=data,robust=0,detail=0,
			         clusters=data[,cluster],
				 theta.des=dependX0)#,Nit=20,step=1.0,notaylor=1)
	pre.theta <- dependX0 %*% twostage.fit$theta
        X <- model.matrix(formula,data)
        Xnam <- colnames(X)
	X <- cbind(X,pre.theta)
	colnames(X)[ncol(X)] <- "pre.theta"
        ww <- fast.reshape(cbind(X,".num"=seq(nrow(X)),".lefttime"=ltimes),varying=c(".num",".lefttime"),id=data[,cluster])
	if (anyNA(ww)) stop("not balanced for predict two stage\n");  
        Prob <- predict.two.stage(twostage.fit,X=ww[,Xnam],
                                  times=ww[,".lefttime1"],times2=ww[,".lefttime2"],
                                  theta=ww$pre.theta)
        P0 <- numeric(nrow(X))
        P0[ww[,".num1"]] <- Prob$St1t2
        P0[ww[,".num2"]] <- Prob$St1t2
        data[,weight.name2] <- P0
	print(summary(P0))
    }
    data[,weight.name] <- pr
    ## }}} 
        
    if (same.cens & missing(cluster)) message("no cluster for same-cens given \n"); 

    if (same.cens & !missing(cluster)) { ## {{{         
###        message("Minimum weights...")
        myord <- order(data[,cluster])
        data <- data[myord,,drop=FALSE]
        id <-  table(data[,cluster])
        if (pairs) {
            gem <- data[,cluster]%in%(names(id)[id==2])
            id <- id[id==2]
            data <- data[gem,]
        }
        d0 <- subset(data,select=c(cluster,weight.name))
        noncens <- with(data,!eval(terms(formula)[[2]][[3]]))
        d0[,"observed."] <- noncens
        timevar <- paste("_",cluster,weight.name,sep="")
        d0[,timevar] <- unlist(lapply(id,seq))
        Wide <- reshape(d0,direction="wide",timevar=timevar,idvar=cluster)
        W <- apply(Wide[,paste(weight.name,1:2,sep=".")],1,
                   function(x) min(x,na.rm=TRUE))
        Wmarg <- d0[,weight.name]
        data[,weight.name] <- 1/Wmarg
        Wmin <- rep(W,id)
        
        obs1only <- rep(with(Wide, observed..1 & (is.na(observed..2) | !observed..2)),id)
        obs2only <- rep(with(Wide, observed..2 & (is.na(observed..1) | !observed..1)),id)
        obsOne <- which(na.omit(obs1only|obs2only))
        obsBoth <- rep(with(Wide, !is.na(observed..1) & !is.na(observed..2) & observed..2 & observed..1),id)

        data[obsBoth,weight.name] <-
            ifelse(noncens[obsBoth],1/Wmin[obsBoth],0)    
        data[obsOne,weight.name] <-
            ifelse(noncens[obsOne],1/Wmarg[obsOne],0)
    } ## }}} 

    if (obs.only)
        data <- data[noncens,,drop=FALSE]
    
    return(data)    
} ## }}} 

##' Internal function.
##' Calculates Inverse Probability of Censoring and Truncation 
##' Weights and adds them to a data.frame
##'
##' @title Inverse Probability of Censoring Weights
##' @param data data frame
##' @param times possible time argument for speciying a maximum value of time tau=max(times), to specify when things are considered censored or not.
##' @param entrytime nam of entry-time for truncation. 
##' @param time name of time variable on data frame.  
##' @param cause name of cause indicator on data frame.
##' @param same.cens For clustered data, should same censoring be assumed and same truncation (bivariate probability calculated as mininum of the marginal probabilities)
##' @param cluster name of clustering variable
##' @param pairs For paired data (e.g. twins) only the complete pairs are returned (With pairs=TRUE)
##' @param strata  name of strata variable to get weights stratified.
##' @param obs.only Return data with uncensored observations only
##' @param cens.formula model for Cox models for truncation and right censoring times.
##' @param cens.code censoring.code 
##' @param pair.cweight Name of weight variable in the new data.frame for right censorig of pairs
##' @param pair.tweight Name of weight variable in the new data.frame for left truncation of pairs
##' @param pair.weight Name of weight variable in the new data.frame for right censoring and left truncation of pairs
##' @param cname Name of weight variable in the new data.frame for right censoring of individuals
##' @param tname Name of weight variable in the new data.frame for left truncation of individuals
##' @param weight.name Name of weight variable in the new data.frame for right censoring and left truncation of individuals
##' @param prec.factor To let tied censoring and truncation times come after the death times. 
##' @param ... Additional arguments to censoring model 
##' @author Thomas Scheike 
##' @examples
##' d <- simnordic.random(3000,delayed=TRUE,ptrunc=0.7,
##'       cordz=0.5,cormz=2,lam0=0.3,country=FALSE)
##' d$strata <- as.numeric(d$country)+(d$zyg=="MZ")*4
##' times <- seq(60,100,by=10)
##' c1 <- comp.risk(Event(time,cause)~1+cluster(id),data=d,cause=1,
##' 	model="fg",times=times,max.clust=NULL,n.sim=0)
##' mm=model.matrix(~-1+zyg,data=d)
##' out1<-random.cif(c1,data=d,cause1=1,cause2=1,same.cens=TRUE,theta.des=mm)
##' summary(out1)
##' pc1 <- predict(c1,X=1,se=0)
##' plot(pc1)
##' 
##' dl <- d[!d$truncated,]
##' dl <- ipw2(dl,cluster="id",same.cens=TRUE,time="time",entrytime="entry",cause="cause",
##'            strata="strata",prec.factor=100)
##' cl <- comp.risk(Event(time,cause)~+1+
##' 		cluster(id),
##' 		data=dl,cause=1,model="fg",
##' 		weights=dl$indi.weights,cens.weights=rep(1,nrow(dl)),
##'             times=times,max.clust=NULL,n.sim=0)
##' pcl <- predict(cl,X=1,se=0)
##' lines(pcl$time,pcl$P1,col=2)
##' mm=model.matrix(~-1+factor(zyg),data=dl)
##' out2<-random.cif(cl,data=dl,cause1=1,cause2=1,theta.des=mm,
##'                  weights=dl$weights,censoring.weights=rep(1,nrow(dl)))
##' summary(out2)
##' @export
ipw2 <- function(data,times=NULL,entrytime=NULL,time="time",cause="cause",
     same.cens=FALSE,cluster=NULL,pairs=FALSE,
     strata=NULL,obs.only=TRUE,cens.formula=NULL,cens.code=0,
     pair.cweight="pcw",pair.tweight="ptw",pair.weight="weights",
     cname="cweights",tname="tweights",weight.name="indi.weights",
     prec.factor=100)
{ ## {{{ 
   ### first calculates weights based on marginal
   ### estimators  of censoring and truncation

## {{{  weights, up at T_i or min(T_i,max(times))
   if (is.null(times)) times <- max(data[,time])
   if (is.null(entrytime)) entry <- rep(0,nrow(data)) else entry <- data[,entrytime]
   mtt <- max(times)
   prec <- .Machine$double.eps * prec.factor
   trunc.model <- cens.model <- NULL ## output of Cox models for entry cens

   if (is.null(cens.formula)) { 
   if (is.null(strata)) { ## {{{ 
	   if (!is.null(entrytime)) {
		   surv.trunc <- survival::survfit(Surv(-data[,time],-entry+prec,rep(1,nrow(data))) ~ 1) 
		   trunc.dist <- summary(surv.trunc)
		   trunc.dist$time <- rev(-trunc.dist$time)
		   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
		   Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),data[,time],strict=TRUE)
		   Lw <- Lfit[,2]
	   } else {Lw <- 1; }
	   if (!is.null(entrytime)) 
	   ud.cens<- survival::survfit(Surv(entry,data[,time],data[,cause]==0)~+1) 
           else ud.cens<- survival::survfit(Surv(data[,time],data[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,data[,time]),strict=TRUE)[,2];
           weights <- 1/(Lw*Gcx); 
	   cweights <-  Gcx; 
	   tweights <-  Lw; 
   ### ## }}} 
   } else { ## {{{ 
	   ### compute for each strata and combine 
	  vstrata <- as.numeric(factor(data[,strata]))
          weights <- rep(1,nrow(data))
          cweights <- rep(1,nrow(data))
          tweights <- rep(1,nrow(data))
	  for (i in unique(vstrata)) { ## {{{ for each strata
	       who <- (vstrata == i)
	       if (sum(who) <= 1) stop(paste("strata",i,"less than 1 observation\n")); 
	   datas <- subset(data,who)
	   entrytimes <- entry[who]
	   if (!is.null(entrytime)) {
		   surv.trunc <- survival::survfit(Surv(-datas[,time],-entrytimes+prec,rep(1,nrow(datas))) ~ +1) 
		   trunc.dist <- summary(surv.trunc)
		   trunc.dist$time <- rev(-trunc.dist$time)
		   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
		   Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),datas[,time],strict=TRUE)
		   Lw <- Lfit[,2]
	   } else {Lw <- 1; }
	   if (!is.null(entrytime)) 
	   ud.cens<- survival::survfit(Surv(entrytimes,datas[,time],datas[,cause]==0)~+1) 
           else ud.cens<- survival::survfit(Surv(entrytimes,datas[,time],datas[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,datas[,time]),strict=TRUE)[,2];
	   weights[who]<-  1/(Lw*Gcx); 
	   cweights[who] <-  Gcx; 
	   tweights[who] <-  Lw; 
          } ## }}} 
   } ## }}} 
   } else { ### cens.formula Cox models  ## {{{
        X <- model.matrix(cens.formula,data=data)[,-1,drop=FALSE]; 

	if (!is.null(entrytime)) {
###		trunc.model <- cox.aalen(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ prop(X)) 
		trunc.model <- survival::coxph(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ X) 
###                baseout <- cens.model$cum 
		baseout <- survival::basehaz(trunc.model,centered=FALSE); 
		baseout <- cbind(rev(-baseout$time),rev(baseout$hazard))
	###
		Lfit <-Cpred(baseout,data[,time],strict=TRUE)[,-1]
		RR<-exp(as.matrix(X) %*% coef(trunc.model))
		Lfit<-exp(-Lfit*RR)
		Lw <- Lfit
        } else {Lw <- 1; }
###
	if (!is.null(entrytime)) 
	cens.model <- survival::coxph(Surv(entrytime,data[,time],data[,cause]==0)~+X) 
        else cens.model <- survival::coxph(Surv(data[,time],data[,cause]==0)~+X) 
        baseout <- survival::basehaz(cens.model,centered=FALSE); 
###        baseout <- cens.model$cum 
	baseout <- cbind(baseout$time,baseout$hazard)
	Gfit<-Cpred(baseout,pmin(mtt,data[,time]),strict=TRUE)[,2];
	RR<-exp(as.matrix(X) %*% coef(cens.model))
	Gfit<-exp(-Gfit*RR)
        weights <- 1/(Lw*Gfit); 
        cweights <- Gfit
        tweights <- Lw
   } ## }}} 

   data[,weight.name] <- weights
   data[,"cw"] <- 1

   if (!is.null(entrytime)) {
   mint <- min(tweights); maxt <- min(tweights) 
   if (mint<0 | mint>1) warning("min(truncation weights) strange, maybe prec.factor should be different\n")
   if (maxt<0 | maxt>1) warning("max(truncation weights) strange, maybe prec.factor should be different\n")
   }

   attr(data,"trunc.model") <- trunc.model
   attr(data,"cens.model") <- cens.model 
   data[,cname] <- cweights
   data[,tname] <- tweights
## }}} 

  observed <- ((data[,time]>mtt & data[,cause]==cens.code)) | (data[,cause]!=cens.code)
  data[,"observed."] <- observed

  if (same.cens & is.null(cluster)) message("no cluster for same-cens given \n"); 
    if (same.cens & !is.null(cluster)) { ## {{{         
###        message("Minimum weights.cens..")
###        message("Maximum weights.trunc..")
        myord <- order(data[,cluster])
        data <- data[myord,,drop=FALSE]
        id <-  table(data[,cluster])
        observed <- ((data[,time]>mtt & data[,cause]==cens.code)) | (data[,cause]!=cens.code)
       d0 <- subset(data,select=c(cluster,cname,tname))
       noncens <- observed
       d0[,"observed."] <-observed
       timevar <- paste("_",cluster,cname,sep="")
       d0[,timevar] <- unlist(lapply(id,seq))
###    print(head(d0))
###    Wide <- reshape(d0,direction="wide",timevar=timevar,idvar=cluster)
       Wide <- fast.reshape(d0,id=cluster)
       ### censoring same cens
       W <- apply(Wide[,paste(cname,1:2,sep="")],1,
                   function(x) min(x)) ## NA when there is just one
###                   function(x) min(x,na.rm=TRUE))
###        Wmarg <- d0[,cname]
###        data[,pair.cweight] <- 0; ##1/Wmarg
###	print(dim(data))
        Wmin <- rep(W,id)
	data[,pair.cweight] <- 1/Wmin

      ### truncation same truncation  
   if (!is.null(entrytime)) { 
        Wt <- apply(Wide[,paste(tname,1:2,sep="")],1,
                   function(x) min(x)) ## NA when there is just one
###                function(x) max(x,na.rm=TRUE))
###        Wtmarg <- d0[,tname]
###        data[,pair.tweight] <- 0; ##1/Wtmarg
        Wtmax <- rep(Wt,id)
	data[,pair.tweight] <- 1/Wtmax
   }

   ## {{{ 
###        obs1only <- rep(with(Wide, observed.1 & (is.na(observed.2) | !observed.2)),id)
###        obs2only <- rep(with(Wide, observed.2 & (is.na(observed.1) | !observed.1)),id)
###        obsOne <- which(na.omit(obs1only|obs2only))
###        obsBoth <- rep(with(Wide, !is.na(observed.1) & !is.na(observed.2) & observed.2 & observed.1),id)
###        
###        data[obsBoth,pair.cweight] <-
###            ifelse(noncens[obsBoth],1/Wmin[obsBoth],0)    
###        data[obsOne,pair.cweight] <-
###            ifelse(noncens[obsOne],1/Wmarg[obsOne],0)
###   if (!is.null(entrytime)) {
###        data[obsBoth,pair.tweight] <-
###            ifelse(noncens[obsBoth],1/Wtmax[obsBoth],0)    
###        data[obsOne,pair.tweight] <-
###            ifelse(noncens[obsOne],1/Wtmarg[obsOne],0)
###   }
   ## }}} 

    if (!is.null(entrytime)) { 
    data[,pair.weight] <- data[,pair.cweight]*data[,pair.tweight]
    } else data[,pair.weight] <- data[,pair.cweight]
    } ## }}} 

    if (obs.only) { data <- data[observed,] } 

    if (pairs) {
          id <-  table(data[,cluster])
          gem <- data[,cluster]%in%(names(id)[id==2])
          id <- id[id==2]
          data <- data[gem,]
    }

   return(data)
} ## }}} 

