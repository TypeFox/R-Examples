
make.formula<-function(names) formula(paste("~",paste(names,collapse="+")))

dimnames.survey.design<-function(x) dimnames(x$variables)
dimnames.svyrep.design<-function(x) dimnames(x$variables)
dimnames.twophase<-function(x) dimnames(x$phase1$sample$variables)

oldsvydesign<-function(ids,probs=NULL,strata=NULL,variables=NULL, fpc=NULL,
                    data=NULL, nest=FALSE, check.strata=!nest,weights=NULL){
 
    .Deprecated("svydesign")
  
    ## less memory-hungry version for sparse tables
    interaction<-function (..., drop = TRUE) {
        args <- list(...)
        narg <- length(args)
        if (narg == 1 && is.list(args[[1]])) {
            args <- args[[1]]
            narg <- length(args)
        }
        
        ls<-sapply(args,function(a) length(levels(a)))
        ans<-do.call("paste",c(lapply(args,as.character),sep="."))
        ans<-factor(ans)
        return(ans)
        
    }


    na.failsafe<-function(object,...){
      if (NCOL(object)==0)
        object
      else na.fail(object)
    }
    
     if(inherits(ids,"formula")) {
	 mf<-substitute(model.frame(ids,data=data,na.action=na.failsafe))   
	 ids<-eval.parent(mf)
	} else if (!is.null(ids))
            ids<-na.fail(data.frame(ids))

     if(inherits(probs,"formula")){
	mf<-substitute(model.frame(probs,data=data,na.action=na.failsafe))
	probs<-eval.parent(mf)
	}
     
     if(inherits(weights,"formula")){
       mf<-substitute(model.frame(weights,data=data,na.action=na.failsafe))
       weights<-eval.parent(mf)
     } else if (!is.null(weights))
         weights<-na.fail(data.frame(weights))
     
     if(!is.null(weights)){
       if (!is.null(probs))
         stop("Can't specify both sampling weights and probabilities")
       else
         probs<-1/weights
     }

      

    if (!is.null(strata)){
      if(inherits(strata,"formula")){
        mf<-substitute(model.frame(strata,data=data, na.action=na.failsafe))
        strata<-eval.parent(mf)
      }
      if(is.list(strata))
        strata<-na.fail(do.call("interaction", strata))
      if (!is.factor(strata))
        strata<-factor(strata)
      has.strata<-TRUE
    } else {
      strata<-factor(rep(1,NROW(ids)))
      has.strata <-FALSE
    }
    
    if (inherits(variables,"formula")){
        mf<-substitute(model.frame(variables,data=data,na.action=na.pass))
        variables <- eval.parent(mf)
    } else if (is.null(variables)){
        variables<-data
    } else
        variables<-data.frame(variables)

    
     if (inherits(fpc,"formula")){
       mf<-substitute(model.frame(fpc,data=data,na.action=na.failsafe))
       fpc<-eval.parent(mf)
       if (length(fpc))
         fpc<-fpc[,1]
     }
     
    if (is.null(ids) || NCOL(ids)==0)
	ids<-data.frame(.id=seq(length=NROW(variables)))

     ## force subclusters nested in clusters
     if (nest && NCOL(ids)>1){
      N<-ncol(ids)
      for(i in 2:(N)){
          ids[,i]<-do.call("interaction", ids[,1:i,drop=TRUE])
      }
    }
     ## force clusters nested in strata
     if (nest && has.strata && NCOL(ids)){
       N<-NCOL(ids)
       for(i in 1:N)
         ids[,i]<-do.call("interaction", list(strata, ids[,i]))
     }

    ## check if clusters nested in strata 
     if (check.strata && nest)
      warning("No point in check.strata=TRUE if nest=TRUE")
    if(check.strata && !is.null(strata) && NCOL(ids)){
       sc<-rowSums(table(ids[,1],strata)>0)
       if(any(sc>1)) stop("Clusters not nested in strata")
    }

    ## Put degrees of freedom (# of PSUs in each stratum) in object, to 
    ## allow subpopulations
    if (NCOL(ids)){
        nPSU<-table(strata[!duplicated(ids[,1])])
    }


    if (!is.null(fpc)){

       if (NCOL(ids)>1){
         if (all(fpc<1))
           warning("FPC is not currently supported for multi-stage sampling")
         else
           stop("Can't compute FPC from population size for multi-stage sampling")
       }
       
       ## Finite population correction: specified per observation
       if (is.numeric(fpc) && length(fpc)==NROW(variables)){
         tbl<-by(fpc,list(strata),unique)
         if (any(sapply(tbl,length)!=1))
           stop("fpc not constant within strata")
         fpc<-data.frame(strata=factor(rownames(tbl),levels=levels(strata)),
                         N=as.vector(tbl))
       }
       ## Now reduced to fpc per stratum
       nstr<-table(strata[!duplicated(ids[[1]])])
       
       if (all(fpc[,2]<=1)){
         fpc[,2]<- nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
       } else if (any(fpc[,2]<nstr[match(as.character(fpc[,1]), names(nstr))]))
         stop("Over 100% sampling in some strata")
       
     }

    ## if FPC specified, but no weights, use it for weights
    if (is.null(probs) && is.null(weights) && !is.null(fpc)){
      pstr<-nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
      probs<-pstr[match(as.character(strata),as.character(fpc[,1]))]
      probs<-as.vector(probs)
    }


    certainty<-rep(FALSE,length(unique(strata)))
    names(certainty)<-as.character(unique(strata))
    if (any(nPSU==1)){
      ## lonely PSUs: are they certainty PSUs?
      if (!is.null(fpc)){
        certainty<- fpc$N < 1.01
        names(certainty)<-as.character(fpc$strata)
      } else if (all(as.vector(probs)<=1)){
        certainty<- !is.na(match(as.character(unique(strata)),as.character(strata)[probs > 0.99]))
        names(certainty)<-as.character(unique(strata))
      } else {
        warning("Some strata have only one PSU and I can't tell if they are certainty PSUs")
      }
      
    } 
    
    if (is.numeric(probs) && length(probs)==1)
        probs<-rep(probs, NROW(variables))
    
    if (length(probs)==0) probs<-rep(1,NROW(variables))
    
    if (NCOL(probs)==1) probs<-data.frame(probs)

    rval<-list(cluster=ids)
    rval$strata<-strata
    rval$has.strata<-has.strata
    rval$prob<- apply(probs,1,prod) 
    rval$allprob<-probs
    rval$call<-match.call()
    rval$variables<-variables
    rval$fpc<-fpc
    rval$certainty<-certainty
    rval$call<-sys.call()
    rval$nPSU<-nPSU
    class(rval)<-"survey.design"
    rval
  }

print.survey.design<-function(x,varnames=FALSE,design.summaries=FALSE,...){
  .svycheck(x)
  n<-NROW(x$cluster)
  if (x$has.strata) cat("Stratified ")
  un<-length(unique(x$cluster[,1]))
  if(n==un){
    cat("Independent Sampling design\n")
    is.independent<-TRUE
  } else {
    cat(NCOL(x$cluster),"- level Cluster Sampling design\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
    is.independent<-FALSE
  }
  print(x$call)
  if (design.summaries){
    cat("Probabilities:\n")
    print(summary(x$prob))
    if(x$has.strata){
      cat("Stratum sizes: \n")
      a<-rbind(obs=table(x$strata),
	       design.PSU=x$nPSU,
               actual.PSU=if(!is.independent || !is.null(x$fpc))
               table(x$strata[!duplicated(x$cluster[,1])]))
      print(a)
    }
    if (!is.null(x$fpc)){
      if (x$has.strata) {
        cat("Population stratum sizes (PSUs): \n")
        print(x$fpc)
      } else {
        cat("Population size (PSUs):",x$fpc[,2],"\n")
      }
    }
  }
  if (varnames){
    cat("Data variables:\n")
    print(names(x$variables))
  }
  invisible(x)
}

"[.survey.design"<-function (x,i, ...){
  
  if (!missing(i)){ 
    if (is.calibrated(x)){
      tmp<-x$prob[i,]
      x$prob<-rep(Inf, length(x$prob))
      x$prob[i,]<-tmp
    } else {
      x$variables<-"[.data.frame"(x$variables,i,...,drop=FALSE)
      x$cluster<-x$cluster[i,,drop=FALSE]
      x$prob<-x$prob[i]
      x$allprob<-x$allprob[i,,drop=FALSE]
      x$strata<-x$strata[i]
    }
  } else {
    x$variables<-x$variables[,...,drop=FALSE]
  }
  
  x
}

"[<-.survey.design"<-function(x, ...,value){
  if (inherits(value, "survey.design"))
    value<-value$variables
  x$variables[...]<-value
  x
}

dim.survey.design<-function(x){
	dim(x$variables)
}

na.fail.survey.design<-function(object,...){
	tmp<-na.fail(object$variables,...)
	object
}

na.omit.survey.design<-function(object,...){
  tmp<-na.omit(object$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object<-object[-omit,]
    object$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

na.exclude.survey.design<-function(object,...){
	tmp<-na.exclude(object$variables,...)
	exclude<-attr(tmp,"na.action")
	if (length(exclude)){
           object<-object[-exclude,]
	   object$variables<-tmp
	   attr(object,"na.action")<-exclude
	}
	object
}


update.survey.design<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$variables[,newnames[j]]<-eval(dots[[j]],object$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}


subset.survey.design<-function(x,subset,...){
        e <- substitute(subset)
        r <- eval(e, x$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
	x$call<-sys.call(-1)
	x
}

summary.survey.design<-function(object,...){
  class(object)<-"summary.survey.design"
  object
}

print.summary.survey.design<-function(x,...){
  y<-x
  class(y)<-"survey.design"
  print(y,varnames=TRUE,design.summaries=TRUE,...)
}	
     
postStratify.survey.design<-function(design, strata, population, partial=FALSE,...){

  if(inherits(strata,"formula")){
    mf<-substitute(model.frame(strata, data=design$variables,na.action=na.fail))
    strata<-eval.parent(mf)
  }
  strata<-as.data.frame(strata)

  sampletable<-xtabs(I(1/design$prob)~.,data=strata)
  sampletable<-as.data.frame(sampletable)

  if (inherits(population,"table"))
    population<-as.data.frame(population)
  else if (is.data.frame(population))
    population$Freq <- as.vector(population$Freq) ##allows Freq computed by tapply()
  else
    stop("population must be a table or dataframe")

  if (!all(names(strata) %in% names(population)))
    stop("Stratifying variables don't match")
  nn<- names(population) %in% names(strata)
  if (sum(!nn)!=1)
    stop("stratifying variables don't match")

  names(population)[which(!nn)]<-"Pop.Freq"
  
  both<-merge(sampletable, population, by=names(strata), all=TRUE)

  samplezero <- both$Freq %in% c(0, NA)
  popzero <- both$Pop.Freq %in% c(0, NA)
  both<-both[!(samplezero & popzero),]
  
  if (any(onlysample<- popzero & !samplezero)){
    print(both[onlysample,])
    stop("Strata in sample absent from population. This Can't Happen")
  }
  if (any(onlypop <- samplezero & !popzero)){
    if (partial){
      both<-both[!onlypop,]
      warning("Some strata absent from sample: ignored")
    } else {
      print(both[onlypop,])
      stop("Some strata absent from sample: use partial=TRUE to ignore them.")
    }
  } 

  reweight<-both$Pop.Freq/both$Freq
  both$label <- do.call("interaction", list(both[,names(strata)]))
  designlabel <- do.call("interaction", strata)
  index<-match(designlabel, both$label)

  attr(index,"oldweights")<-1/design$prob
  design$prob<-design$prob/reweight[index]
  attr(index,"weights")<-1/design$prob
  design$postStrata<-c(design$postStrata,list(index))
  
  ## Do we need to iterate here a la raking to get design strata
  ## and post-strata both balanced?
  design$call<-sys.call(-1)
  
  design
}


svyCprod<-function(x, strata, psu, fpc, nPSU, certainty=NULL, postStrata=NULL,
                   lonely.psu=getOption("survey.lonely.psu")){

  x<-as.matrix(x)
  n<-NROW(x)

  ## Remove post-stratum means, which may cut across PSUs
  if(!is.null(postStrata)){
    for (psvar in postStrata){
      if (inherits(psvar, "greg_calibration") || inherits(psvar, "raking"))
        stop("rake() and calibrate() not supported for old-style design objects")
      psw<-attr(psvar,"weights")
      psmeans<-rowsum(x/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
      x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
    }
  }

  ##First collapse over PSUs

  if (is.null(strata)){
    strata<-rep("1",n)
    if (!is.null(nPSU))
        names(nPSU)<-"1"
  }
  else
    strata<-as.character(strata) ##can't use factors as indices in for()'

  if (is.null(certainty)){
    certainty<-rep(FALSE,length(strata))
    names(certainty)<-strata
  }
  
  if (!is.null(psu)){
    x<-rowsum(x, psu, reorder=FALSE)
    strata<-strata[!duplicated(psu)]
    n<-NROW(x)
  }
  
  if (!is.null(nPSU)){
      obsn<-table(strata)
      dropped<-nPSU[match(names(obsn),names(nPSU))]-obsn
      if(sum(dropped)){
        xtra<-matrix(0,ncol=NCOL(x),nrow=sum(dropped))
        strata<-c(strata,rep(names(dropped),dropped))
      	if(is.matrix(x))
	   x<-rbind(x,xtra)
        else
	   x<-c(x,xtra)
        n<-NROW(x)
      }
  } else obsn<-table(strata)

  if(is.null(strata)){
      x<-t(t(x)-colMeans(x))
  } else {
      strata.means<-drop(rowsum(x,strata, reorder=FALSE))/drop(rowsum(rep(1,n),strata, reorder=FALSE))
      if (!is.matrix(strata.means))
          strata.means<-matrix(strata.means, ncol=NCOL(x))
      x<- x- strata.means[ match(strata, unique(strata)),,drop=FALSE]
  }
  
  p<-NCOL(x)
  v<-matrix(0,p,p)
  
  ss<-unique(strata)
  for(s in ss){
      this.stratum <- strata %in% s
      
      ## original number of PSUs in this stratum 
      ## before missing data/subsetting
      this.n <-nPSU[match(s,names(nPSU))]
      
      this.df <- this.n/(this.n-1)	
      
      if (is.null(fpc))
          this.fpc <- 1
      else{
          this.fpc <- fpc[,2][ fpc[,1]==as.character(s)]
          this.fpc <- (this.fpc - this.n)/this.fpc
      }
      
      xs<-x[this.stratum,,drop=FALSE]

      this.certain<-certainty[names(certainty) %in% s]
      
      ## stratum with only 1 design cluster leads to undefined variance
      lonely.psu<-match.arg(lonely.psu, c("remove","adjust","fail",
                                          "certainty","average"))
      if (this.n==1 && !this.certain){
        this.df<-1
        if (lonely.psu=="fail")
          stop("Stratum ",s, " has only one sampling unit.")
        else if (lonely.psu!="certainty")
          warning("Stratum ",s, " has only one sampling unit.")
        if (lonely.psu=="adjust")
          xs<-strata.means[match(s,ss),,drop=FALSE]
      } else if (obsn[match(s,names(obsn))]==1 && !this.certain){
        ## stratum with only 1 cluster left after subsetting 
        warning("Stratum ",s," has only one PSU in this subset.")
        if (lonely.psu=="adjust")
          xs<-strata.means[match(s,ss),,drop=FALSE]
      }
      ## add it up
      if (!this.certain)
        v<-v+crossprod(xs)*this.df*this.fpc
    }
  if (lonely.psu=="average"){
    v<- v/(1-mean(obsn==1 & !certainty))
  }
  v
}



svymean<-function(x, design,na.rm=FALSE,...){
  .svycheck(design)
  UseMethod("svymean",design)
}

svymean.survey.design<-function(x,design, na.rm=FALSE,deff=FALSE,...){

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  }
  else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$variables)
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
            design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-svyCprod(x*pweights/psum,design$strata,design$cluster[[1]], design$fpc,
              design$nPSU,design$certainty, design$postStrata)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
    nobs<-NROW(design$cluster)
    vsrs<-svyvar(x,design,na.rm=na.rm)/nobs
    vsrs<-vsrs*(psum-nobs)/psum
    attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}


print.svystat<-function(x,...){
    vv<-attr(x,"var")
    if (is.matrix(vv))
        m<-cbind(x,sqrt(diag(vv)))
    else
        m<-cbind(x,sqrt(vv))
    hasdeff<-!is.null(attr(x,"deff"))
    if (hasdeff) {
        m<-cbind(m,deff(x))
        colnames(m)<-c(attr(x,"statistic"),"SE","DEff")
    } else {
        colnames(m)<-c(attr(x,"statistic"),"SE")
    }
    printCoefmat(m)
}

as.data.frame.svystat<-function(x,...){
  rval<-data.frame(statistic=coef(x),SE=SE(x))
  names(rval)[1]<-attr(x,"statistic")
  if (!is.null(attr(x,"deff")))
    rval<-cbind(rval,deff=deff(x))
  rval
}

coef.svystat<-function(object,...){
  attr(object,"statistic")<-NULL
  attr(object,"deff")<-NULL
  attr(object,"var")<-NULL
  unclass(object)
}

vcov.svystat<-function(object,...){
  as.matrix(attr(object,"var"))
}

SE.svystat<-function(object,...){
 v<-vcov(object)
 if (!is.matrix(v) || NCOL(v)==1) sqrt(v) else sqrt(diag(v))
}

deff <- function(object,quietly=FALSE,...) UseMethod("deff")

deff.default <- function(object, quietly=FALSE,...){
  rval<-attr(object,"deff")
  if (is.null(rval)) { 
    if(!quietly)
      warning("object has no design effect information")
  } else rval<-diag(as.matrix(rval))
  rval
}

cv<-function(object,...) UseMethod("cv")

cv.default<-function(object, warn=TRUE, ...){
  rval<-SE(object)/coef(object)
  if (warn && any(coef(object)<0,na.rm=TRUE)) warning("CV may not be useful for negative statistics")
  rval
}


svytotal<-function(x,design,na.rm=FALSE,...){
  .svycheck(design)
  UseMethod("svytotal",design)
}
svytotal.survey.design<-function(x,design, na.rm=FALSE, deff=FALSE,...){

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
    if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  } else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, design$variables)

  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }

  N<-sum(1/design$prob)
  m <- svymean(x, design, na.rm=na.rm)
  total<-m*N
  attr(total, "var")<-v<-svyCprod(x/design$prob,design$strata,
                                  design$cluster[[1]], design$fpc,
                                  design$nPSU,design$certainty,design$postStrata)
  attr(total,"statistic")<-"total"
  if (is.character(deff) || deff){
    vsrs<-svyvar(x,design)*sum(weights(design)^2)
    vsrs<-vsrs*(N-NROW(design$cluster))/N
    attr(total,"deff")<-v/vsrs
  }
  return(total)
}

svyvar<-function(x, design, na.rm=FALSE,...){
  .svycheck(design)
  UseMethod("svyvar",design)
}
svyvar.survey.design<-function(x, design, na.rm=FALSE,...){
    
	if (inherits(x,"formula"))
            x<-model.frame(x,model.frame(design),na.action=na.pass)
	else if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        
	n<-sum(weights(design,"sampling")!=0)
	xbar<-svymean(x,design, na.rm=na.rm)
	if(NCOL(x)==1) {
            x<-x-xbar
            v<-svymean(x*x*n/(n-1),design, na.rm=na.rm)
            attr(v,"statistic")<-"variance"
            return(v)
	}
	x<-t(t(x)-xbar)
	p<-NCOL(x)
	a<-matrix(rep(x,p),ncol=p*p)
	b<-x[,rep(1:p,each=p)]
        ## Kish uses the n-1 divisor, so it affects design effects
	v<-svymean(a*b*n/(n-1),design, na.rm=na.rm)
	vv<-matrix(v,ncol=p)
        dimnames(vv)<-list(names(xbar),names(xbar))
        attr(vv,"var")<-attr(v,"var")
        attr(vv,"statistic")<-"variance"
        class(vv)<-c("svyvar","svystat")
        vv
    }

print.svyvar<-function (x,  covariance=FALSE, ...) 
{
    if(!is.matrix(x)) NextMethod()
    
    vv <- attr(x, "var")
    if (covariance){
      nms<-outer(rownames(x),colnames(x),paste,sep=":")
      m<-cbind(as.vector(x), sqrt(diag(vv)))
      rownames(m)<-nms
    } else{
      ii <- which(diag(sqrt(length(x)))>0)
      m <- cbind(x[ii], sqrt(diag(vv))[ii])
    }
    colnames(m) <- c(attr(x, "statistic"), "SE")
    printCoefmat(m)
}

as.matrix.svyvar<-function(x,...) unclass(x)

svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)

svyquantile.survey.design<-function(x,design,quantiles,alpha=0.05,
                                    ci=FALSE, method="linear",f=1,
                                    interval.type=c("Wald","score","betaWald"),
                                    na.rm=FALSE,se=ci, ties=c("discrete","rounded"), df=Inf,...){
    if (inherits(x,"formula"))
      x<-model.frame(x ,model.frame(design), na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, model.frame(design,na.action=na.pass))
    
    if (na.rm){
        nas<-rowSums(is.na(x))
        design<-design[nas==0,]
        if (length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
        else
          x[nas>0,]<-0
      }
   

    w<-weights(design)

    if (is.null(df)){
      qcrit<-function(p, lower.tail=TRUE) qt(p, df=degf(design), lower.tail=lower.tail)
    } else if(df==Inf){
      qcrit <- function(p,lower.tail=TRUE) qnorm(p,lower.tail=lower.tail)
    } else {
      qcrit <- function(p,lower.tail=TRUE) qt(p,df=df,lower.tail=lower.tail)
    }

    
    computeQuantiles<-function(xx,p=quantiles){
      if (any(is.na(x))) return(NA*p)
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method=method,f=f,
                     yleft=min(xx),yright=max(xx),ties=min) 
      cdf(p)
    }
    
    computeQuantilesRounded<-function(xx,p=quantiles){
      if (any(is.na(xx))) return(NA*p)
      ww<-rowsum(w,xx,reorder=TRUE)
      xx<-sort(unique(xx))
      cum.w <- cumsum(ww)/sum(ww)
      cdf <- approxfun(cum.w, xx, method = method, f = f, 
                       yleft = min(xx), yright = max(xx),ties=min)
      cdf(p)
    }
      
    
    
    computeScoreCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
   
      U<-function(theta){ ((xx>theta)-(1-p))}
        
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }
      
      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),
                qlimit=qcrit(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),
                qlimit=qcrit(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }
    
    computePCI<-function(se,alpha,p){
      if (interval.type=="Wald"){
        p.up<-p+qcrit(alpha/2,lower.tail=FALSE)*se
        p.low<-p+qcrit(alpha/2,lower.tail=TRUE)*se
        c(p.low,p.up)
      } else if (interval.type=="betaWald"){
        n.eff <- (p*(1-p))/(se^2)
        n.eff <- n.eff * ( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
        p.up<-qbeta(1-alpha/2, n.eff*p+1, n.eff*(1-p))
        p.low<-qbeta(alpha/2,  n.eff*p, n.eff*(1-p)+1)
        c(p.low,p.up)
      }
      
    }
    
    computeWaldCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      theta0<-computeQuantiles(xx,p)
      U<- ((xx>theta0)-(1-p))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      approx(cum.w,xx[oo],xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx),ties=min)$y 
      
    }
    
    computeWaldCIRounded<-function(xx,p){
      if(any(is.na(xx))) return(c(NA,NA))
        theta0<-computeQuantilesRounded(xx,p)
        U<- ((xx>theta0)-(1-p))
        ww<-rowsum(w,xx, reorder=TRUE)
        uxx <- sort(unique(xx))
        wtest<-svymean(U,design)
        p.ci<-computePCI(SE(wtest),alpha,p)
        p.low<-p.ci[1]
        p.up<-p.ci[2]
        oo<-order(xx)
        cum.w<-cumsum(ww)/sum(ww)
        approx(cum.w,uxx,xout=c(p.low,p.up), method=method,f=f,
               yleft=min(xx),yright=max(xx),ties=min)$y 
        
      }

    ties<-match.arg(ties)
    computeQ<-switch(ties, discrete=computeQuantiles,rounded=computeQuantilesRounded)
    
    if (!is.null(dim(x)))
        rval<-t(matrix(apply(x,2,computeQ),nrow=length(quantiles),
                       dimnames=list(as.character(round(quantiles,2)),colnames(x))))
    else
      rval<-computeQ(x)
    
    if (!ci & !se) return(rval)
    
    interval.type<-match.arg(interval.type)
    
    computeCI<-switch(paste(interval.type,ties,sep="."), score.discrete=computeScoreCI,
                            score.rounded=stop("ties=\"rounded\" not available with interval.type=\"score\""),
                            Wald.rounded=computeWaldCIRounded,
                            betaWald.rounded=computeWaldCIRounded,
                            Wald.discrete=computeWaldCI,
                            betaWald.discrete=computeWaldCI)
    
    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
                   as.character(round(quantiles,2)),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    if (ci)
      rval<-list(quantiles=rval,CIs=cis)
    else
      rval<-list(quantiles=rval)
    
    if (is.null(dim(x)))
        ses<-(cis[2,]-cis[1,])/(2*qcrit(alpha/2,lower.tail=FALSE))
    else
        ses<-(cis[2,,]-cis[1,,])/(2*qcrit(alpha/2,lower.tail=FALSE))
    attr(rval,"SE")<-ses
    class(rval)<-"svyquantile"
    rval
  }

SE.svyquantile<-function(object,...){
    attr(object,"SE")
}

vcov.svyquantile<-function(object,...){
  se<-SE(object)
  if (is.null(se)) stop("no uncertainty information present")
  v<-matrix(NA,length(se),length(se))
  warning("Only diagonal of vcov() available")
  diag(v)<-se
  v
}

coef.svyquantile<-function(object,...){
  rval<-as.vector(object$quantiles)
  if(ncol(object$quantiles)==1)
    names(rval)<-rownames(object$quantiles)
  else if (nrow(object$quantiles)==1)
    names(rval)<-colnames(object$quantiles)
  else names(rval)<-t(outer(colnames(object$quantiles),
                            rownames(object$quantiles),
                            paste,sep=":"))
  rval
}

print.svyquantile<-function(x,...){
    print(list(quantiles=x$quantiles, CIs=x$CIs))
}

coef.svyratio<-function(object,...,drop=TRUE){
  if (!drop) return(object$ratio)
  cf<-as.vector(object$ratio)
  nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
  names(cf)<-nms
  cf
}
 
SE.svyratio<-function(object,...,drop=TRUE){
  if(!drop) return(sqrt(object$var))
  se<-as.vector(sqrt(object$var))
  nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
  names(se)<-nms
  se
}

svyratio<-function(numerator,denominator, design,...){
  .svycheck(design)
  UseMethod("svyratio",design)
}

svyratio.survey.design<-function(numerator, denominator, design,...){

    if (inherits(numerator,"formula"))
		numerator<-model.frame(numerator,design$variables)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
		denominator<-model.frame(denominator,design$variables)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    allstats<-svytotal(all,design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-svyCprod(r*1/design$prob, design$strata, design$cluster[[1]], design$fpc,
                            design$nPSU, design$certainty,design$postStrata)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    rval$call<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }

print.svyratio_separate<-function(x,...){
  cat("Stratified ratio estimate: ")
  if (!is.null(x$call))
    print(x$call)
  else if (!is.null(attr(x,"call")))
    print(attr(x$call))
  for(r in x$ratios) {
    print(r)
  }
  invisible(x)
}

print.svyratio<-function(x,...){
  cat("Ratio estimator: ")
  if (!is.null(x$call))
    print(x$call)
  else if(!is.null(attr(x,"call")))
    print(attr(x,"call"))
  cat("Ratios=\n")
  print(x$ratio)
  cat("SEs=\n")
  print(sqrt(x$var))
  invisible(NULL)
}

predict.svyratio<-function(object, total, se=TRUE,...){
  if (se)
    return(list(total=object$ratio*total,se=sqrt(object$var)*total))
  else
    return(object$ratio*total)
}

predict.svyratio_separate<-function(object, total, se=TRUE,...){

  if (length(total)!=length(object$ratios))
    stop("Number of strata differ in ratio object and totals.")
  if (!is.null(names(total)) && !is.null(levels(object$strata))){
    if (!setequal(names(total), levels(object$strata)))
      warning("Names of strata differ in ratio object and totals")
    else if (!all(names(total)==levels(object$strata))){
      warning("Reordering supplied totals to make their names match the ratio object")
      total<-total[match(names(total),levels(object$strata))]
    }
  }
  totals<-mapply(predict, object=object$ratios, total=total,se=se,...,SIMPLIFY=FALSE)

  if(se){
    rval<-totals[[1]]$total
    v<-totals[[1]]$se^2
    for(ti in totals[-1]) {
      rval<-rval+ti$total
      v<-v+ti$se^2
    }
    list(total=rval,se=sqrt(v))
  } else {
    rval<-totals[[1]]
    for (ti in totals[-1]) rval<-rval+ti
    rval
  }

}


cv.svyratio<-function(object,...){
  sqrt(object$var)/object$ratio
}

svytable<-function(formula, design, ...){
    UseMethod("svytable",design)
}

svytable.survey.design<-function(formula, design, Ntotal=NULL, round=FALSE,...){
  
  if (!inherits(design,"survey.design")) stop("design must be a survey design")
  weights<-1/design$prob
  
  ## unstratified or unadjusted
  if (length(Ntotal)<=1 || !design$has.strata){
    if (length(formula)==3)
      tblcall<-bquote(xtabs(I(weights*.(formula[[2]]))~.(formula[[3]]), data=model.frame(design),...))
       else
         tblcall<-bquote(xtabs(weights~.(formula[[2]]), data=model.frame(design),...))
    tbl<-eval(tblcall)
    if (!is.null(Ntotal)) {
      if(length(formula)==3)
        tbl<-tbl/sum(Ntotal)
      else
        tbl<-tbl*sum(Ntotal)/sum(tbl)
    }
    if (round)
      tbl<-round(tbl)
    attr(tbl,"call")<-match.call()
    class(tbl)<-c("svytable",class(tbl))
    return(tbl)
  }
  ## adjusted and stratified
  if (length(formula)==3)
    tblcall<-bquote(xtabs(I(weights*.(formula[[2]]))~design$strata[,1]+.(formula[[3]]), data=model.frame(design),...))
  else
    tblcall<-bquote(xtabs(weights~design$strata[,1]+.(formula[[2]]), data=model.frame(design),...))
  
  tbl<-eval(tblcall)
  
  ss<-match(sort(unique(design$strata[,1])), Ntotal[,1])
  dm<-dim(tbl)
  layer<-prod(dm[-1])
  tbl<-sweep(tbl,1,Ntotal[ss, 2]/apply(tbl,1,sum),"*")
  tbl<-apply(tbl, 2:length(dm), sum)
  if (round)
    tbl<-round(tbl)
  class(tbl)<-c("svytable","xtabs", "table")
  attr(tbl, "call")<-match.call()
  tbl
}

svycoxph<-function(formula,design,subset=NULL,...){
  .svycheck(design)
  UseMethod("svycoxph",design)
}

svycoxph.survey.design<-function(formula,design,subset=NULL,...){
    subset<-substitute(subset)
    subset<-eval(subset, model.frame(design),parent.frame())
    if (!is.null(subset))
        design<-design[subset,]

    if(any(weights(design)<0)) stop("weights must be non-negative")
    
    require(survival) || stop("Needs the survival package")
    data<-model.frame(design)
    
    g<-match.call()
    g$formula<-eval.parent(g$formula)
    g$design<-NULL
    g$var<-NULL
    if (is.null(g$weights))
      g$weights<-quote(.survey.prob.weights)
    else
      g$weights<-bquote(.survey.prob.weights*.(g$weights))
    g[[1]]<-quote(coxph)
    g$data<-quote(data)
    g$subset<-quote(.survey.prob.weights>0)
    g$model <- TRUE
    
    ##need to rescale weights for stability
    data$.survey.prob.weights<-(1/design$prob)/mean(1/design$prob)
    if (!all(all.vars(formula) %in% names(data))) 
        stop("all variables must be in design= argument")
    g<-with(list(data=data), eval(g))
    g$call<-match.call()
    g$call[[1]]<-as.name(.Generic)
    g$printcall<-sys.call(-1)
    g$printcall[[1]]<-as.name(.Generic)
    class(g)<-c("svycoxph", class(g))
    g$survey.design<-design
    
    nas<-g$na.action
    if (length(nas))
        design<-design[-nas,]

    dbeta.subset<-resid(g,"dfbeta",weighted=TRUE)
    if (nrow(design)==NROW(dbeta.subset)){
      dbeta<-as.matrix(dbeta.subset)
    } else {
      dbeta<-matrix(0,ncol=NCOL(dbeta.subset),nrow=nrow(design))
      dbeta[is.finite(design$prob),]<-dbeta.subset
    }
    g$inv.info<-g$var
    
    if (inherits(design,"survey.design2"))
      g$var<-svyrecvar(dbeta, design$cluster,
                    design$strata, design$fpc,
                    postStrata=design$postStrata)
    else if (inherits(design, "twophase"))
      g$var<-twophasevar(dbeta, design)
    else if(inherits(design, "twophase2"))
      g$var<-twophase2var(dbeta, design)
    else if(inherits(design, "pps"))
      g$var<-ppsvar(dbeta,design)
    else
      g$var<-svyCprod(dbeta, design$strata,
                      design$cluster[[1]], design$fpc,design$nPSU,
                      design$certainty,design$postStrata)
    
    g$wald.test<-coef(g)%*%solve(g$var,coef(g))
    g$ll<-g$loglik
    g$loglik<-NULL
    g$rscore<-NULL
    g$score<-NA
    g$degf.resid<-degf(design)-length(coef(g)[!is.na(coef(g))])+1
    
    g
}


model.frame.svycoxph<-function(formula,...){
    f<-formula$call
    env <- environment(formula(formula))
    if (is.null(env)) 
        env <- parent.frame()
    f[[1]]<-as.name("model.frame")
    f$data<-quote(data)
    f$design<-NULL
    f$method<-f$control<-f$singular.ok<-f$model<-f$x<-f$y<-f$iter<-NULL
    f$formula<-formula(formula)
    if (is.null(f$weights))
        f$weights<-quote(.survey.prob.weights)
    else 
        f$weights<-bquote(.survey.prob.weights*.(f$weights))
    design<-formula$survey.design
    data<-model.frame(design)
    data$.survey.prob.weights<-(1/design$prob)/sum(1/design$prob)
    with(list(data=data), eval(f))
}

model.matrix.svycoxph<-function (object, data = NULL, contrast.arg = object$contrasts, 
    ...) 
{
    if (!is.null(object[["x"]])) 
        object[["x"]]
    else {
        if (is.null(data)) 
            data <- model.frame(object, ...)
        else data <- model.frame(object, data = data, ...)
        Terms <- object$terms
        attr(Terms, "intercept") <- 1
        strats <- attr(Terms, "specials")$strata
        cluster <- attr(Terms, "specials")$cluster
        dropx <- NULL
        if (length(cluster)) {
            tempc <- untangle.specials(Terms, "cluster", 1:10)
            ord <- attr(Terms, "order")[tempc$terms]
            if (any(ord > 1)) 
                stop("Cluster can not be used in an interaction")
            dropx <- tempc$terms
        }
        if (length(strats)) {
            temp <- untangle.specials(Terms, "strata", 1)
            dropx <- c(dropx, temp$terms)
        }
        if (length(dropx)) {
            newTerms <- Terms[-dropx]
            X <- model.matrix(newTerms, data, contrasts = contrast.arg)
        }
        else {
            newTerms <- Terms
            X <- model.matrix(Terms, data, contrasts = contrast.arg)
        }
        X
    }
}

print.svycoxph<-function(x,...){
    print(x$survey.design, varnames=FALSE, design.summaries=FALSE,...)
##    x$call<-x$printcall
    NextMethod()
}

summary.svycoxph<-function(object,...){
    print(object$survey.design,varnames=FALSE, design.summaries=FALSE,...)
##    object$call<-object$printcall
    NextMethod()
}

survfit.svycoxph<-function(object,...){
    stop("No survfit method for survey models")
}
extractAIC.svycoxph<-function(fit,...){
    stop("No AIC for survey models")
}

anova.svycoxph<-function(object,...){
    stop("No anova method for survey models")
}

svyglm<-function(formula, design, ...){
  .svycheck(design)
  UseMethod("svyglm",design)
}

svyglm.survey.design<-function(formula,design,subset=NULL,...){

      subset<-substitute(subset)
      subset<-eval(subset, model.frame(design), parent.frame())
      if (!is.null(subset))
        design<-design[subset,]
      
      data<-model.frame(design)

      g<-match.call()
      g$formula<-eval.parent(g$formula)
      g$design<-NULL
      g$var<-NULL
      if (is.null(g$weights))
        g$weights<-quote(.survey.prob.weights)
      else 
        g$weights<-bquote(.survey.prob.weights*.(g$weights))
      g$data<-quote(data)
      g[[1]]<-quote(glm)      

      ##need to rescale weights for stability in binomial
      data$.survey.prob.weights<-(1/design$prob)/mean(1/design$prob)
      if (!all(all.vars(formula) %in% names(data))) 
	stop("all variables must be in design= argument")
      g<-with(list(data=data), eval(g))
      g$naive.cov<-summary(g)$cov.unscaled
      
      nas<-g$na.action
      if (length(nas))
	design<-design[-nas,]

      g$cov.unscaled<-svy.varcoef(g,design)
      g$df.residual <- degf(design)+1-length(coef(g)[!is.na(coef(g))])
      
      class(g)<-c("svyglm",class(g))
      g$call<-sys.call()
      g$call[[1]]<-as.name(.Generic)
      if(!("formula" %in% names(g$call))) {
        if (is.null(names(g$call)))
          i<-1
        else
          i<-min(which(names(g$call)[-1]==""))
        names(g$call)[i+1]<-"formula"
      }
      g$survey.design<-design 
      g
}

print.svyglm<-function(x,...){
  print(x$survey.design, varnames=FALSE, design.summaries=FALSE,...)
  NextMethod()

}

coef.svyglm<-function(object,...,na.rm=TRUE) {
  beta<-object$coefficients
  if (!na.rm || length(beta)==object$rank)
    beta
  else
    beta[object$qr$pivot[1:object$rank]]
}

vcov.svyglm<-function(object,...) {
  v<-object$cov.unscaled
  dimnames(v)<-list(names(coef(object)),names(coef(object)))
  v
}


svy.varcoef<-function(glm.object,design){
    Ainv<-summary(glm.object)$cov.unscaled
    estfun<-model.matrix(glm.object)*resid(glm.object,"working")*glm.object$weights
    if (glm.object$rank<NCOL(estfun)){
      estfun<-estfun[,glm.object$qr$pivot[1:glm.object$rank]]
    }
    naa<-glm.object$na.action
    ## the design may still have rows with weight zero for missing values
    ## if there are weights or calibration. model.matrix will have removed them
    if (length(naa) && (NROW(estfun)!=nrow(design) )){
      if ((length(naa)+NROW(estfun))!=nrow(design) )
        stop("length mismatch: this can't happen.")
      n<-nrow(design)     
      inx <- (1:n)[-naa]
      ee <- matrix(0,nrow=n,ncol=NCOL(estfun))
      ee[inx,]<-estfun
      estfun<-ee
    }

    if (inherits(design,"survey.design2"))
      svyrecvar(estfun%*%Ainv,design$cluster,design$strata,design$fpc,postStrata=design$postStrata)
    else if (inherits(design, "twophase"))
      twophasevar(estfun%*%Ainv, design)
    else if (inherits(design, "twophase2"))
      twophase2var(estfun%*%Ainv, design)
    else if (inherits(design, "pps"))
      ppsvar(estfun%*%Ainv, design)
    else
      svyCprod(estfun%*%Ainv,design$strata,design$cluster[[1]],design$fpc, design$nPSU,
                  design$certainty,design$postStrata)
  }

residuals.svyglm<-function(object,type = c("deviance", "pearson", "working", 
    "response", "partial"),...){
	type<-match.arg(type)
	if (type=="pearson"){
   	   y <- object$y
	   mu <- object$fitted.values
    	   wts <- object$prior.weights
           pwts<- 1/object$survey.design$prob
           pwts<- pwts/mean(pwts)
           ## missing values in calibrated/post-stratified designs
           ## the rows will still be in the design object but not in the model
           if (length(naa<-object$na.action) && (length(pwts)!=length(wts))){
             if(length(naa)+length(wts) != length(pwts))
               stop("length mismatch: this can't happen.")
             inx<-(1:length(pwts))[-naa]
           } else inx<-1:length(pwts)
           r<-numeric(length(pwts))
	   r[inx]<-(y - mu) * sqrt(wts/pwts[inx])/(sqrt(object$family$variance(mu)))
	   if (is.null(object$na.action)) 
        	r
    	   else 
	        naresid(object$na.action, r)
	} else 
		NextMethod()

}

summary.svyglm<-function (object, correlation = FALSE, df.resid=NULL,...) 
{
    Qr <- object$qr
    est.disp <- TRUE
    if (is.null(df.resid))
      df.r <- object$df.residual
    else
      df.r<-df.resid
    
    dispersion<-svyvar(resid(object,"pearson"), object$survey.design,
                       na.rm=TRUE)
    
    coef.p <- coef(object)
    covmat<-vcov(object)
    dimnames(covmat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
    }
    else if (df.r > 0) {
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
            "Pr(>|t|)"))
    }
    else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        aic = object$aic, coefficients = coef.table, dispersion = dispersion, 
        df = c(object$rank, df.r,NCOL(Qr$qr)), cov.unscaled = covmat, 
        cov.scaled = covmat))
    if (correlation) {
        dd <- sqrt(diag(covmat))
        ans$correlation <- covmat/outer(dd, dd)
    }
    ans$aliased<-is.na(coef(object,na.rm=FALSE))
    ans$survey.design<-list(call=object$survey.design$call)
    class(ans) <- c("summary.svyglm","summary.glm")
    return(ans)
}


logLik.svyglm<-function(object,...){
   warning("svyglm not fitted by maximum likelihood.")
   object$deviance
}

AIC.svyglm<-function(object,...,k=2){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),extractAIC,k=k))
    } else {
	   extractAIC(object,k=k)
    }
}
extractAIC.svyglm<-function(fit,scale,k=2,...){
	if (length(attr(terms(fit),"factors"))){
	    r<-regTermTest(fit, delete.response(formula(fit)), method="LRT")
	    deltabar<-mean(r$lambda)
	} else {
	    r<-list(lambda=0)
	    deltabar<-NaN
	}
	d<-fit$deviance
	c(eff.p=sum(r$lambda), AIC=d+k*sum(r$lambda),deltabar=deltabar)
}

extractAIC.svrepglm<-extractAIC.svyglm

BIC.svyglm<-function(object,...,maximal){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),dBIC,modelM=maximal))
    } else {
	   dBIC(object,modelM=maximal)
    }
	
	}
	
dBIC<-function(modela,modelM){
	pm<-modela$rank
	pM<-modelM$rank	

	if (any(!(names(coef(modela))%in% names(coef(modelM))))){
		stop("coefficients in model but not in maximal model")
		}
	index<-!(names(coef(modelM))%in% names(coef(modela)))
	n<-1+modela$df.null	
	if(any(index)){
		wald<-coef(modelM)[index]%*%solve(vcov(modelM)[index,index],coef(modelM)[index])
		detDelta<-det(solve(modelM$naive.cov[index,index,drop=FALSE],modelM$cov.unscaled[index,index,drop=FALSE]))
		dbar<-detDelta^(1/(pM-pm))
		nstar<-n/dbar	
	}else {
		wald<-0
		detDelta<-1
		dbar<-1
		nstar=NaN
		}
	c(p=pm, BIC=wald+pm*log(n)+log(detDelta)+deviance(modelM),neff=nstar)
	}	


confint.svyglm<-function(object,parm,level=0.95,method=c("Wald","likelihood"),ddf=Inf,...){
  method<-match.arg(method)
  if(method=="Wald"){
    tlevel<-pt(qnorm(level),df=ddf)
    return(confint.default(object,parm=parm,level=tlevel,...))
  }
  pnames <- names(coef(object))
  if (missing(parm)) 
    parm <- seq_along(pnames)
  else if (is.character(parm))
    parm <- match(parm, pnames, nomatch = 0)
  lambda<-diag(object$cov.unscaled[parm,parm,drop=FALSE])/diag(object$naive.cov[parm,parm,drop=FALSE])
  if(is.null(ddf)) ddf<-object$df.residual
  if (ddf==Inf)
    level<- 1-2*pnorm(qnorm((1-level)/2)*sqrt(lambda))
  else {
    level<- 1-2*pnorm(qt((1-level)/2,df=ddf)*sqrt(lambda))
  }
  rval<-vector("list",length(parm))
  for(i in 1:length(parm)){
    rval[[i]]<-NextMethod(object=object,parm=parm[i],level=level[i],...)
  }
  names(rval)<-pnames[parm]
  if (length(rval)==1)
    rval<-rval[[1]]
  else
    rval<-do.call(rbind,rval)
  attr(rval,"levels")<-level
  rval
}


svymle<-function(loglike, gradient=NULL, design, formulas,
                 start=NULL, control=list(maxit=1000),
                 na.action="na.fail", method=NULL,...){
  if(is.null(method))
    method<-if(is.null(gradient)) "Nelder-Mead" else "nlm"
  
  if (!inherits(design,"survey.design")) 
	stop("design is not a survey.design")
  weights<-weights(design)
  wtotal<-sum(weights)
 
  if (is.null(control$fnscale))
      control$fnscale <- -wtotal/length(weights)
  if (inherits(design, "twophase"))
    data<-design$phase1$sample$variables
  else 
    data<-design$variables

## Get the response variable
  nms<-names(formulas)
  if (nms[1]==""){
	if (inherits(formulas[[1]],"formula"))
	  y<-eval.parent(model.frame(formulas[[1]],data=data,na.action=na.pass))
	else
	  y<-eval(y,data,parent.frame())
	formulas[1]<-NULL
	if (FALSE && NCOL(y)>1) stop("Y has more than one column")
    }   else {
  	## one formula must have response
	has.response<-sapply(formulas,length)==3
	if (sum(has.response)!=1) stop("Need a response variable")
	ff<-formulas[[which(has.response)]]
	ff[[3]]<-1
	y<-eval.parent(model.frame(ff,data=data,na.action=na.pass))
	formulas[[which(has.response)]]<-delete.response(terms(formulas[[which(has.response)]]))
        nms<-c("",nms)
  }

  if(length(which(nms==""))>1) stop("Formulas must have names")

  mf<-vector("list",length(formulas))
  vnms <- unique(do.call(c, lapply(formulas, all.vars)))
  uformula <- make.formula(vnms)
  mf <- model.frame(uformula, data=data,na.action=na.pass)
  mf <- cbind(`(Response)`=y, mf)
  mf<-mf[,!duplicated(colnames(mf)),drop=FALSE]
  
  mf<-get(na.action)(mf)  
  nas<-attr(mf,"na.action")
  if (length(nas))
    design<-design[-nas,]
  weights<-1/design$prob
  wtotal<-sum(weights)
  
  Y<-mf[,1]
  mm<-lapply(formulas,model.matrix, data=mf)

  ## parameter names
  parnms<-lapply(mm,colnames)
  for(i in 1:length(parnms))
	parnms[[i]]<-paste(nms[i+1],parnms[[i]],sep=".")
  parnms<-unlist(parnms)

  # maps position in theta to model matrices
  np<-c(0,cumsum(sapply(mm,NCOL)))


  objectivefn<-function(theta,...){
     args<-vector("list",length(nms))
     args[[1]]<-Y
     for(i in 2:length(nms))
	args[[i]]<-mm[[i-1]]%*%theta[(np[i-1]+1):np[i]]
     names(args)<-nms
     args<-c(args, ...)
     sum(do.call("loglike",args)*weights)
  }

  if (is.null(gradient)) {
     grad<-NULL
  } else {  
     fnargs<-names(formals(loglike))[-1]
     grargs<-names(formals(gradient))[-1]
     if(!identical(fnargs,grargs))
       stop("loglike and gradient have different arguments.")
     reorder<-na.omit(match(grargs,nms[-1]))
     grad<-function(theta,...){
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       rval<-NULL
       tmp<-do.call("gradient",args)
       for(i in reorder){
	   rval<-c(rval, colSums(as.matrix(tmp[,i]*weights*mm[[i]])))
	}
       drop(rval)
     }
  }

  theta0<-numeric(np[length(np)])
  if (is.list(start))
      st<-do.call("c",start)
  else
      st<-start

  if (length(st)==length(theta0)) {
	theta0<-st
  } else {
	stop("starting values wrong length")
  }

  if (method=="nlm"){
      ff<-function(theta){
          rval<- -objectivefn(theta)
          if (is.na(rval)) rval<- -Inf
          attr(rval,"gradient")<- -grad(theta)
          rval
      }
      rval<-nlm(ff, theta0,hessian=TRUE)
      if (rval$code>3) warning("nlm did not converge")
      rval$par<-rval$estimate
  } else {
      rval<-optim(theta0, objectivefn, grad,control=control,
              hessian=TRUE,method=method,...)
      if (rval$conv!=0) warning("optim did not converge")
  }

  


  names(rval$par)<-parnms
  dimnames(rval$hessian)<-list(parnms,parnms)

  if (is.null(gradient)) {
	rval$invinf<-solve(-rval$hessian)
	rval$scores<-NULL
	rval$sandwich<-NULL
    }  else {
       theta<-rval$par
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       deta<-do.call("gradient",args)
       rval$scores<-NULL
       for(i in reorder)
       	 rval$scores<-cbind(rval$scores,deta[,i]*weights*mm[[i]])

       rval$invinf<-solve(-rval$hessian)
       dimnames(rval$invinf)<-list(parnms,parnms)

       db<-rval$scores%*%rval$invinf
       if (inherits(design,"survey.design2"))
         rval$sandwich<-svyrecvar(db,design$cluster,design$strata, design$fpc, 
                               postStrata=design$postStrata)
       else if (inherits(design, "twophase"))
         rval$sandwich<-twophasevar(db,design)
       else
         rval$sandwich<-svyCprod(db,design$strata,design$cluster[[1]],
                                 design$fpc, design$nPSU,
                                 design$certainty, design$postStrata)
       dimnames(rval$sandwich)<-list(parnms,parnms)
     }
  rval$call<-match.call()
  rval$design<-design
  class(rval)<-"svymle"
  rval

}


svymleOLD<-function(loglike, gradient=NULL, design, formulas,
                 start=NULL, control=list(maxit=1000),
                 na.action="na.fail", method=NULL,...){
  if(is.null(method))
    method<-if(is.null(gradient)) "Nelder-Mead" else "nlm"
  
  if (!inherits(design,"survey.design")) 
	stop("design is not a survey.design")

  weights<-1/design$prob
  wtotal<-sum(weights)
  if (is.null(control$fnscale))
      control$fnscale<- -wtotal
  if (inherits(design, "twophase"))
    data<-design$phase1$sample$variables
  else 
    data<-design$variables

## Get the response variable
  nms<-names(formulas)
  if (nms[1]==""){
	if (inherits(formulas[[1]],"formula"))
	  y<-eval.parent(model.frame(formulas[[1]],data=data,na.action=na.pass))
	else
	  y<-eval(y,data,parent.frame())
	formulas[1]<-NULL
	if (FALSE && NCOL(y)>1) stop("Y has more than one column")
    }   else {
  	## one formula must have response
	has.response<-sapply(formulas,length)==3
	if (sum(has.response)!=1) stop("Need a response variable")
	ff<-formulas[[which(has.response)]]
	ff[[3]]<-1
	y<-eval.parent(model.frame(ff,data=data,na.action=na.pass))
	formulas[[which(has.response)]]<-delete.response(terms(formulas[[which(has.response)]]))
        nms<-c("",nms)
  }

  if(length(which(nms==""))>1) stop("Formulas must have names")
  
  
  mf<-vector("list",length(formulas))
  for(i in 1:length(formulas)){
	mf[[i]]<-eval.parent(model.frame(formulas[[i]], data=data, na.action=na.pass))
	}
  notnulls<-sapply(mf,function(mfi) NCOL(mfi)!=0)
  mf<-mf[notnulls]
  if (any(notnulls))
    mf<-as.data.frame(do.call("cbind",c(y,mf)))
  else
    mf<-y
  names(mf)[1]<-"(Response)"
  mf<-mf[,!duplicated(colnames(mf)),drop=FALSE]

  mf<-get(na.action)(mf)  
  nas<-attr(mf,"na.action")
  if (length(nas))
	design<-design[-nas,]

  Y<-mf[,1]
  mm<-lapply(formulas,model.matrix, data=mf)

  ## parameter names
  parnms<-lapply(mm,colnames)
  for(i in 1:length(parnms))
	parnms[[i]]<-paste(nms[i+1],parnms[[i]],sep=".")
  parnms<-unlist(parnms)

  # maps position in theta to model matrices
  np<-c(0,cumsum(sapply(mm,NCOL)))


  objectivefn<-function(theta,...){
     args<-vector("list",length(nms))
     args[[1]]<-Y
     for(i in 2:length(nms))
	args[[i]]<-mm[[i-1]]%*%theta[(np[i-1]+1):np[i]]
     names(args)<-nms
     args<-c(args, ...)
     sum(do.call("loglike",args)*weights)
  }

  if (is.null(gradient)) {
     grad<-NULL
  } else {  
     fnargs<-names(formals(loglike))[-1]
     grargs<-names(formals(gradient))[-1]
     if(!identical(fnargs,grargs))
       stop("loglike and gradient have different arguments.")
     reorder<-na.omit(match(grargs,nms[-1]))
     grad<-function(theta,...){
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       rval<-NULL
       tmp<-do.call("gradient",args)
       for(i in reorder){
	   rval<-c(rval, colSums(as.matrix(tmp[,i]*weights*mm[[i]])))
	}
       drop(rval)
     }
  }

  theta0<-numeric(np[length(np)])
  if (is.list(start))
      st<-do.call("c",start)
  else
      st<-start

  if (length(st)==length(theta0)) {
	theta0<-st
  } else {
	stop("starting values wrong length")
  }

  if (method=="nlm"){
      ff<-function(theta){
          rval<- -objectivefn(theta)
          if (is.na(rval)) rval<- -Inf
          attr(rval,"grad")<- -grad(theta)
          rval
      }
      rval<-nlm(ff, theta0,hessian=TRUE)
      if (rval$code>3) warning("nlm did not converge")
      rval$par<-rval$estimate
  } else {
      rval<-optim(theta0, objectivefn, grad,control=control,
              hessian=TRUE,method=method,...)
      if (rval$conv!=0) warning("optim did not converge")
  }

  


  names(rval$par)<-parnms
  dimnames(rval$hessian)<-list(parnms,parnms)

  if (is.null(gradient)) {
	rval$invinf<-solve(-rval$hessian)
	rval$scores<-NULL
	rval$sandwich<-NULL
    }  else {
       theta<-rval$par
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       deta<-do.call("gradient",args)
       rval$scores<-NULL
       for(i in reorder)
       	 rval$scores<-cbind(rval$scores,deta[,i]*weights*mm[[i]])

       rval$invinf<-solve(-rval$hessian)
       dimnames(rval$invinf)<-list(parnms,parnms)

       db<-rval$scores%*%rval$invinf
       if (inherits(design,"survey.design2"))
         rval$sandwich<-svyrecvar(db,design$cluster,design$strata, design$fpc, 
                               postStrata=design$postStrata)
       else if (inherits(design, "twophase"))
         rval$sandwich<-twophasevar(db,design)
       else
         rval$sandwich<-svyCprod(db,design$strata,design$cluster[[1]],
                                 design$fpc, design$nPSU,
                                 design$certainty, design$postStrata)
       dimnames(rval$sandwich)<-list(parnms,parnms)
     }
  rval$call<-match.call()
  rval$design<-design
  class(rval)<-"svymle"
  rval

}

coef.svymle<-function(object,...) object$par

vcov.svymle<-function(object,stderr=c("robust","model"),...) {
    stderr<-match.arg(stderr)
    if (stderr=="robust"){
	rval<-object$sandwich
	if (is.null(rval)) {
		p<-length(coef(object))
		rval<-matrix(NA,p,p)
	}
    } else {
        rval<-object$invinf*mean(1/object$design$prob)
    }
    rval
}


print.svymle<-function(x,...){
  cat("Survey-sampled mle: \n")
  print(x$call)
  cat("Coef:  \n")
  print(x$par)
}

summary.svymle<-function(object,stderr=c("robust","model"),...){
    cat("Survey-sampled mle: \n")
    print(object$call)
    stderr<-match.arg(stderr)
    tbl<-data.frame(Coef=coef(object),SE=sqrt(diag(vcov(object,stderr=stderr))))
    tbl$p.value<-format.pval(2*(1-pnorm(abs(tbl$Coef/tbl$SE))), digits=3,eps=0.001)
    print(tbl)
    print(object$design)
}

model.frame.survey.design<-function(formula,...,drop=TRUE){
  formula$variables
}
model.frame.svyrep.design<-function(formula,...){
  formula$variables
}
model.frame.survey.design2<-function(formula,...){
  formula$variables
}

.onLoad<-function(...){
  if (is.null(getOption("survey.lonely.psu")))
    options(survey.lonely.psu="fail")
  if (is.null(getOption("survey.ultimate.cluster")))
    options(survey.ultimate.cluster=FALSE)
  if (is.null(getOption("survey.want.obsolete")))
    options(survey.want.obsolete=FALSE)
  if (is.null(getOption("survey.adjust.domain.lonely")))
    options(survey.adjust.domain.lonely=FALSE)
  if (is.null(getOption("survey.drop.replicates")))
      options(survey.drop.replicates=TRUE)
  if (is.null(getOption("survey.multicore")))
    options(survey.multicore=FALSE)
  if (is.null(getOption("survey.replicates.mse")))
    options(survey.replicates.mse=FALSE)
}


predterms<-function(object,se=FALSE,terms=NULL){
  tt<-terms(object)
  n <- length(object$residuals)
  p <- object$rank
  p1 <- seq_len(p)
  piv <- object$qr$pivot[p1]
  beta<-coef(object)
  X<-mm<-model.matrix(object)
  aa <- attr(mm, "assign")
  ll <- attr(tt, "term.labels")
  hasintercept <- attr(tt, "intercept") > 0L
  if (hasintercept) 
    ll <- c("(Intercept)", ll)
  aaa <- factor(aa, labels = ll)
  asgn <- split(order(aa), aaa)
  if (hasintercept) {
    asgn$"(Intercept)" <- NULL
    }
  avx <- colMeans(mm)
  termsconst <- sum(avx[piv] * beta[piv])
  nterms <- length(asgn)
  ip <- matrix(ncol = nterms, nrow = NROW(X))
  if (nterms > 0) {
    predictor <- matrix(ncol = nterms, nrow = NROW(X))
    dimnames(predictor) <- list(rownames(X), names(asgn))
    
    if (hasintercept) 
      X <- sweep(X, 2L, avx, check.margin = FALSE)
    unpiv <- rep.int(0L, NCOL(X))
    unpiv[piv] <- p1
    for (i in seq.int(1L, nterms, length.out = nterms)) {
      iipiv <- asgn[[i]]
      ii <- unpiv[iipiv]
      iipiv[ii == 0L] <- 0L
      predictor[, i] <- if (any(iipiv > 0L)) 
        X[, iipiv, drop = FALSE] %*% beta[iipiv]
      else 0
      if (se){
        ip[,i]<-if (any(iipiv > 0L)) 
          rowSums(as.matrix(X[, iipiv, drop = FALSE] %*% vcov(object)[ii,ii]) * X[, iipiv, drop = FALSE]) else 0
      }
    }
    if (!is.null(terms)) {
      predictor <- predictor[, terms, drop = FALSE]
      if (se) 
        ip <- ip[, terms, drop = FALSE]
    }
  }
  else {
    predictor <- ip <- matrix(0, n, 0)
  }
  attr(predictor, "constant") <- if (hasintercept) 
    termsconst
  else 0
  if(se)
          dimnames(ip)<-dimnames(predictor)
  if (se) list(fit=predictor,se.fit=sqrt(ip)) else predictor
}


predict.svyglm <- function(object, newdata=NULL, total=NULL,
                           type = c("link", "response","terms"),
                           se.fit=(type!="terms"),
                           vcov=FALSE,...){
    if(is.null(newdata))
      newdata<-model.frame(object$survey.design)
    type<-match.arg(type)
    if (type=="terms")
      return(predterms(object,se=se.fit,...))
    tt<-delete.response(terms(formula(object)))
    mf<-model.frame(tt,data=newdata)
    mm<-model.matrix(tt,mf)
    if (!is.null(total) && attr(tt,"intercept")){
        mm[,attr(tt,"intercept")]<-mm[,attr(tt,"intercept")]*total
    }
    eta<-drop(mm %*% coef(object))
    d<-drop(object$family$mu.eta(eta))
    eta<-switch(type, link=eta, response=object$family$linkinv(eta))
    if(se.fit){
        if(vcov){
            vv<-mm %*% vcov(object) %*% t(mm)
            attr(eta,"var")<-switch(type,
                                    link=vv,
                                    response=d*(t(vv*d)))
        } else {
            ## FIXME make this more efficient
            vv<-drop(rowSums((mm %*% vcov(object)) * mm))
            attr(eta,"var")<-switch(type,
                                    link=vv,
                                    response=drop(d*(t(vv*d))))
        }
    }
    attr(eta,"statistic")<-type
    class(eta)<-"svystat"
    eta
    }
    
