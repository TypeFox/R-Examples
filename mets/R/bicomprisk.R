##' Estimation of concordance in bivariate competing risks data
##'
##' @title Estimation of concordance in bivariate competing risks data
##' @param formula Formula with left-hand-side being a \code{Event} object (see example below) and the left-hand-side specying the covariate structure
##' @param data Data frame
##' @param cause Causes (default (1,1)) for which to estimate the bivariate cumulative incidence
##' @param cens The censoring code
##' @param causes causes
##' @param indiv indiv
##' @param strata Strata
##' @param id Clustering variable
##' @param num num
##' @param max.clust max number of clusters in comp.risk call for iid decompostion, max.clust=NULL uses all clusters otherwise rougher grouping.
##' @param marg marginal cumulative incidence to make stanard errors for same clusters for subsequent use in casewise.test()
##' @param se.clusters to specify clusters for standard errors. Either a vector of cluster indices or a column name in \code{data}. Defaults to the \code{id} variable.
##' @param wname name of additonal weight used for paired competing risks data. 
##' @param prodlim prodlim to use prodlim estimator (Aalen-Johansen) rather than IPCW weighted estimator based on comp.risk function.These are equivalent in the case of no covariates. These esimators are the same in the case of stratified fitting. 
##' @param messages Control amount of output
##' @param model Type of competing risk model (default is Fine-Gray model "fg", see comp.risk). 
##' @param return.data Should data be returned (skipping modeling)
##' @param uniform to compute uniform standard errors for concordance estimates based on resampling.
##' @param conservative for conservative standard errors, recommended for larger data-sets.
##' @param resample.iid to return iid residual processes for further computations such as tests. 
##' @param ... Additional arguments to comp.risk function 
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' 
##' ## Simulated data example 
##' prt <- simnordic.random(2000,delayed=TRUE,ptrunc=0.7,
##'		      cordz=0.5,cormz=2,lam0=0.3)
##' ## Bivariate competing risk, concordance estimates
##' p11 <- bicomprisk(Event(time,cause)~strata(zyg)+id(id),
##'                   data=prt,cause=c(1,1))
##'  
##' p11mz <- p11$model$"MZ"
##' ## Concordance
##' plot(p11mz,ylim=c(0,0.1)); 
##' 
##' ## entry time, truncation weighting 
##' prtl <- prt[!prt$truncated,]
##' prt2 <- ipw2(prtl,cluster="id",same.cens=TRUE,
##'      time="time",entrytime="entry",cause="cause",
##'      strata="zyg",pairs=TRUE,obs.only=FALSE)
##' prt2$ppt <- prt2$ptw
##'
##' p11t <- bicomprisk(Event(time,cause)~strata(zyg)+id(id),
##' 		  data=prt2,cause=c(1,1),wname="ppt",times=50:90)
##' ## Concordance
##' p11tmz <- p11t$model$"MZ"
##' p11tdz <- p11t$model$"DZ"
##' lines(p11tmz$time,p11tmz$P1,col=2)
##' 
##' ### other weighting procedure 
##' prt2 <- ipw2(prt,cluster="id",same.cens=TRUE,
##'      time="time",cause="cause",
##'      pairs=TRUE,strata="zyg",obs.only=TRUE)
##'
##' prt22 <- fast.reshape(prt2,id="id")
##'
##' prt22$event <- (prt22$cause1==1)*(prt22$cause2==1)*1
##' prt22$timel <- pmax(prt22$time1,prt22$time2)
##' ipwc <- comp.risk(Event(timel,event)~-1+factor(zyg1),
##'   data=prt22,cause=1,n.sim=0,model="rcif2",times=50:90,
##'   weights=prt22$weights1,cens.weights=rep(1,nrow(prt22)))
##'
##' p11wmz <- ipwc$cum[,2] 
##' p11wdz <- ipwc$cum[,3]
##' lines(ipwc$cum[,1],p11wmz,col=3)
##' 
bicomprisk <- function(formula, data, cause=c(1,1), cens=0, causes, indiv, 
 strata=NULL, id,num, max.clust=1000, marg=NULL,se.clusters=NULL,wname=NULL,
 prodlim=FALSE,messages=TRUE,model,return.data=0,uniform=0,conservative=1,resample.iid=1,...) {
## {{{ 

  mycall <- match.call()
  formulaId <- unlist(Specials(formula,"id"))
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- unlist(Specials(formula,"strata"))
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)

  ### ts  11/10 
  ## if (substr(as.character(formula)[2],1,4)=="Hist") {
  ##     stop("Since version : The left hand side of the formula must be specified as 
  ##     Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
  ## }


  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  } 
  if (missing(id)) stop("Missing 'id' variable")
  
  ### setting up cluster level for iid decomposition
  if (is.null(se.clusters) & is.null(marg))  lse.clusters <- data[,c(id)]
  else { 
    if (is.null(se.clusters)) {
      lse.clusters <- marg$clusters
###      if (!is.null(max.clust)) {   }
    } else {
      if (is.character(se.clusters)) {
        lse.clusters <- data[,se.clusters]
      } else {
        lse.clusters <- se.clusters
      }
    }
  }
  if (length(lse.clusters)!=nrow(data)) stop("'se.clusters' and 'data' does not match!")
  se.clusters.call <- se.clusters

  data <- data.frame(cbind(data,lse.clusters))

  timevar <- terms(formula)[[2]]
  ##  hh <- with(data,eval(timevar))
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  
  if (!is.null(strata)) {
    fac <- interaction(data[,strata])
    dd <- split(data,fac)
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq_len(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        idx <- which(fac==names(dd)[i])
        mycall$se.clusters <- lse.clusters[idx]
        mycall$formula <- formula
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("bicomprisk.strata","twinlm.strata")
      res$N <- length(dd)
      return(res)
    }
  }

  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  ### adds weights, 
  if (!is.null(wname)) covars <- c(covars,wname)

  indiv2 <- covars2 <- NULL 
  
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  ##which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  oldreshape <- 0
  if (oldreshape==1) sep="." else sep=""
  timevar2 <- paste(timevar,1:2,sep=sep)
  causes2 <- paste(causes,1:2,sep=sep)
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=sep)
  for (i in seq_len(length(indiv)))
  indiv2 <- c(indiv2, paste(indiv[i],1:2,sep=sep))

  if (oldreshape==1)
  ww0 <- reshape(data[,c(timevar,causes,covars,indiv,id,num,"lse.clusters")],
         direction="wide",idvar=id,timevar=num)[,c(id,"lse.clusters.1",timevar2,causes2,indiv2,covars2)]
  else 
  ww0 <- fast.reshape(data[,c(timevar,causes,covars,indiv,id,num,"lse.clusters")],id=id,num=data$num,labelnum=TRUE)[,c(id,"lse.clusters1",timevar2,causes2,indiv2,covars2)]
  ww0 <- na.omit(ww0)
 
  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[1]]
  
  ## {{{ (i,j) causes 
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2])
  if (length(idx2)>0) {
      status[idx2] <- 1
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,max)
  }

  ##(0,0), (0,j)
  idx2 <- which(ww0[,causes2[1]]==cens & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(ic,0), (ic,j)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[1]]
  }

  ##(i,0)
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens)
  if (length(idx2)>0) {
      status[idx2] <- 0
      time[idx2] <- ww0[idx2,timevar2[2]]
  }

  ##(ic,jc)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- apply(ww0[idx2,timevar2[1:2],drop=FALSE],1,min)
  }

  ##(0,jc),(i,jc)
  idx2 <- which((ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[1]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  if (length(idx2)>0) {
      status[idx2] <- 2
      time[idx2] <- ww0[idx2,timevar2[2]]
  }
  
  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2])
  names(mydata) <- c(timevar,causes,covars,indiv2)
  ## }}}

  if (return.data==2) return(list(data=mydata)) else {
  if (!prodlim) {
    ff <- paste("Event(",timevar,",",causes,",cens.code=",cens,") ~ 1",sep="")
    if (!is.null(wname)) covars <- covars[-which(covars %in% c(wname))]
    if (length(c(covars,indiv))>0) {
     xx <- c(covars,indiv2)
     for (i in seq_len(length(xx))) 
	     xx[i] <- paste("const(",xx[i],")",sep="")
      ff <- paste(c(ff,xx),collapse="+")
      if (missing(model)) model <- "fg"      
    }
    if (missing(model)) model <- "fg"
    ### clusters for iid construction
    lse.clusters <- NULL
    if (!is.null(se.clusters.call)) {
        lse.clusters <- ww0[,"lse.clusters1"]
    }

    if (!(is.null(wname))) mydata <- ipw2(mydata,time=timevar,cause=causes) #,cens.code=cens)

    if (is.null(wname)) {
	    add<-comp.risk(as.formula(ff),data=mydata,
	    cause=1,n.sim=0,resample.iid=resample.iid,model=model,conservative=conservative,
	    clusters=lse.clusters, max.clust=max.clust,...)
    } else { 
	    add<-comp.risk(as.formula(ff),data=mydata,
	    cause=1,n.sim=0,resample.iid=resample.iid,model=model,conservative=conservative,
	    clusters=lse.clusters, max.clust=max.clust,
	    weights=mydata[,wname]*mydata$indi.weights,cens.weights=rep(1,nrow(mydata)),...)
    }

    padd <- predict(add,X=1,se=1,uniform=uniform,resample.iid=resample.iid)
    padd$cluster.names <- lse.clusters
  } else {
      if (!requireNamespace("prodlim",quietly=TRUE)) stop("prodlim requested but not installed")
    ff <- as.formula(paste("Hist(",timevar,",",causes,")~",paste(c("1",covars,indiv2),collapse="+")))
    padd <- prodlim::prodlim(ff,data=mydata)
  }
###  class(padd) <- c("bicomprisk",class(padd))
 if (return.data==1) return(list(comp.risk=padd,data=mydata)) else return(padd)  
  }
## }}} 
}

## plot.bicomprisk <- function(x,add=FALSE,...) {
##   if ("predict.timereg"%in%class(a)) {    
##     if (!add) { plot.predict.timereg(x,...) }
##     else {
##       with(x,lines(time,P1,...))
##     }    
##   } else {
##     plot(x,...)
##   }
##   return(invisible(x))
##}
