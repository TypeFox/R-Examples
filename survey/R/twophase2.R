##
## Constructing cov(R_i,R_j)/pi^*_ij, or \check{\check{\Delta}}_ij in Sarndal's notation
## We use this form because it can be sparse and because it is easy to combine
## multistage and multiphase sampling.
##
Dcheck_strat<-function(strata, prob){
    strata<-as.numeric(strata) ## for ave()
    n<-length(strata)
    rval<-matrix(0, n,n)
    sampsize<-ave(strata,strata,FUN=length)
    strats<-unique(strata)
    for(strat in strats){
        these <- strata == strat
        rval[these,these]<- -(1-prob[these])/(sampsize[these]-1)
    }
    diag(rval)<-(1-prob)
    rval
}

Dcheck_multi<-function(id,strata,probs){
   nstage<-NCOL(id)
   rval<-matrix(0,NROW(id),NROW(id))
   for(stage in 1:nstage){
       uid<-!duplicated(id[,stage])
       idx<-match(id[,stage],id[uid,stage])
       this_stage<-Dcheck_strat(strata[uid,stage],probs[uid,stage])[idx,idx]
       rval<- twophaseDcheck(rval, this_stage)
     }
   rval
 }

Dcheck_subset<-function(strata, subset, prob, withreplacement){
    strata<-as.numeric(strata) ## for ave()
    N<-length(strata)
    n<-sum(subset)
    rval<-matrix(0, n,n)
    sampsize<-ave(strata,strata,FUN=length)
    strats<-unique(strata)
    if (!withreplacement){
      for(strat in strats){
        these <- strata == strat
        ithese<-which(these)
        rval[these[subset],these[subset]]<- -(1-prob[ithese[subset]])/(sampsize[ithese[subset]]-1)
      }
    }
    diag(rval)<-(1-prob[subset])
    rval
  }

Dcheck_multi_subset<-function(id,strata,subset,probs,withreplacement){
   nstage<-NCOL(id)
   n<-sum(subset)
   rval<-matrix(0,n,n)
   if (all(probs==1) && withreplacement)
     return(as(diag(n),"sparseMatrix"))
   for(stage in 1:nstage){
       uid<-!duplicated(id[,stage])
       insubset<-rowsum(as.integer(subset),id[,stage],reorder=FALSE)>0
       idx<-match(id[subset,stage],id[subset,stage][uid])
       this_stage<-Dcheck_subset(strata[uid,stage],insubset,probs[uid,stage],withreplacement)[idx,idx]
       rval<- twophaseDcheck(rval, this_stage)
     }
   rval
 }
twophaseDcheck<-function(Dcheck1,Dcheck2){
  as(-Dcheck1*Dcheck2+Dcheck1+Dcheck2,"sparseMatrix")
}

make_covmat<-function(design1,design2,subset){
  require("Matrix",quietly=TRUE) || stop("These designs require the Matrix package")
  withreplacement<-is.null(design1$fpc$popsize)
  phase1<-Dcheck_multi_subset(design1$cluster, design1$strata, subset, design1$allprob, withreplacement)
  phase2<-Dcheck_multi(design2$cluster, design2$strata, design2$allprob)
  dcheck<-twophaseDcheck(phase1,phase2)
  list(phase1=phase1,phase2=phase2,full=dcheck)
}

##
## Based on twophase(), so it computes some stuff that is no longer necessary.
## Will be pruned in the future.
##
twophase2<-function(id,strata=NULL, probs=NULL, fpc=NULL,
                   subset, data){

  d1<-svydesign(ids=id[[1]],strata=strata[[1]],weights=NULL,
                probs=probs[[1]],fpc=fpc[[1]],data=data)

  if(inherits(subset,"formula"))
    subset<-eval.parent(model.frame(subset,data=data,na.action=na.pass))[[1]]

  if(!is.logical(subset) && sort(unique(subset))==c(0,1))
      subset<-as.logical(subset)

  if (any(is.na(subset))) stop("missing values in 'subset'")
  
  d1s<-svydesign(ids=id[[1]],strata=strata[[1]],weights=NULL,
                probs=probs[[1]],fpc=fpc[[1]],data=data[subset,])
  d1s$prob<-d1$prob[subset]
  d1s$allprob<-d1$allprob[subset,,drop=FALSE]
  
  ## work out phase-two fpc
  if (is.null(fpc[[2]])){
    complete.vars<-names(data)[apply(data, 2, function(v) all(!is.na(v)))]
    if (all(c(all.vars(id[[2]]), all.vars(strata[[2]])) %in% complete.vars)){
      dfpc<-svydesign(ids=id[[2]], strata=strata[[2]], data=data, probs=NULL)
      popsize<-mapply(function(s,i) ave(!duplicated(i),s,FUN=sum), dfpc$strata, dfpc$cluster)
      rm(dfpc)
    } else {
      warning("Second-stage fpc not specified and not computable")
      popsize<-NULL
    }
  } else popsize<-NULL

  d2<-svydesign(ids=id[[2]], strata=strata[[2]], probs=probs[[2]],
                weights=NULL, fpc=fpc[[2]], data=data[subset,])

  ## ugly hack to get nicer labels
  if(!is.null(fpc[[2]])){
    d2call<-bquote(svydesign(ids=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                               fpc=.(fpc[[2]])))
  } else{
    d2call<-bquote(svydesign(ids=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                               fpc=`*phase1*`))
  }
  for(i in names(d2call)[-1])
    d2call[[i]]<-d2call[[i]]
  d2$call<-d2call
  d1call<-bquote(svydesign(ids=.(id[[1]]), strata=.(strata[[1]]), probs=.(probs[[1]]),
                               fpc=.(fpc[[1]])))
  for(i in names(d1call)[-1])
    d1call[[i]]<-d1call[[i]]
  d1$call<-d1call


  ## Add phase 2 fpc and probs if they were computed rather than specified.
  if (!is.null(popsize))
    d2$fpc<-as.fpc(popsize[subset,,drop=FALSE],d2$strata,d2$cluster)
  if(is.null(probs[[2]]) && !is.null(d2$fpc$popsize)){
    d2$allprob<-1/weights(d2$fpc,final=FALSE)
    d2$prob<-apply(as.data.frame(d2$allprob),1,prod)
  }
  
  d2$variables<-NULL
  deltacheck<-make_covmat(d1,d2, subset)

  rval<-list(phase1=list(full=d1,sample=d1s),
             phase2=d2,
             subset=subset, dcheck=deltacheck)
  rval$prob<-rval$phase1$sample$prob

  ## Are phase 2 PSUs the same as Phase 1 USUs, or smaller?
  rval$samescale<- !any(duplicated(d1s$cluster[,NCOL(d1s$cluster)][!duplicated(d2$cluster[,1])]))
  
  ## For each phase 1 sampling unit, need probability of being represented
  ## at phase 2.
  nunique<-function(x) sum(!duplicated(x))
  m<-NCOL(rval$phase1$sample$cluster)
  if(d2$has.strata){
      if (inherits(strata[[2]],"formula"))
          sa<-eval(attr(terms(strata[[2]]),"variables")[[2]],d1$variables)
      else
          sa<-d1$strata[,1]
      cm<-rval$phase1$full$cluster[,m]
      if (nunique(sa)!=nunique(sa[subset]))
          stop("Some phase-2 strata have zero sampling fraction")
      rval$usu<-ave(cm[subset],sa[subset],FUN=nunique)/ave(cm,sa,FUN=nunique)[subset]
  } else {
      rval$usu<-drop(with(rval$phase1$sample,ave(cluster[,m], strata[,m], FUN=nunique))/rval$phase1$full$fpc$sampsize[rval$subset])
  }

  if (length(rval$phase1$sample$prob)==length(d2$prob))
    rval$prob<-rval$phase1$sample$prob*d2$prob
  else{
    rval$prob<-rep(Inf,length(rval$phase1$sample$prob))
    rval$prob[subset]<-rval$prob[subset]*d2$prob
  }
  rval$call<-sys.call()
  class(rval) <- c("twophase2","survey.design")
  rval
}

print.twophase2<-function(x,...){
  cat("Two-phase sparse-matrix design:\n ")
  print(x$call)
  cat("Phase 1:\n")
  print(x$phase1$full)
  cat("Phase 2:\n")
  print(x$phase2)
  invisible(x)
}

summary.twophase2<-function(object,...){
  class(object)<-"summary.twophase2"
  object
}

print.summary.twophase2<-function(x,...,varnames=TRUE){
  cat("Two-phase sparse-matrix design:\n ")
  print(x$call)
   cat("Phase 1:\n")
  print(x$phase1$full,design.summaries=TRUE,varnames=FALSE)
  cat("Phase 2:\n")
  print(x$phase2,design.summaries=TRUE, varnames=FALSE)
  if (varnames){
    cat("Data variables:\n")
    print(names(x$phase1$full$variables))
  }
  invisible(x)
}

twophase2var<-function(x,design){
  ## calibration is allowed at phase one or phase two,
  ## but not for clusters within a phase
  postStrata2<-design$phase2$postStrata
  postStrata1<-design$phase1$full$postStrata
  if (is.null(postStrata1) && is.null(postStrata2)){
    rval<-htvar.matrix(x,design$dcheck$full)
    ph2<-htvar.matrix(x,design$dcheck$phase2)
    attr(rval,"phases")<-list(phase1=rval-ph2,phase2=ph2)
    return(rval)
  }
  if (!is.null(postStrata1)){
    ##phase 1 calibration
    ## x is size of phase-2 sample,need to expand to allow calibration.
    y<-matrix(0,ncol=ncol(x),nrow=length(design$subset))
    y[design$subset,]<-x
    for (psvar in postStrata1){
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage==0){
          ## G-calibration at population level
          y<-qr.resid(psvar$qr,y/psvar$w)*psvar$w
        } else {
          ## G-calibration within clusters
          stop("calibration within clusters not allowed for two-phase designs")
        }
      } else {
        ## ordinary post-stratification
        psw<-attr(psvar, "weights")
        postStrata<-as.factor(psvar)
        psmeans<-rowsum(y/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
        y<- y-psmeans[match(psvar,sort(unique(psvar))),]*psw
        x1<-y[design$subset,,drop=FALSE]
      }
    }
  } else x1<-x
  phase1var<-htvar.matrix(x1,design$dcheck$full)-htvar.matrix(x1,design$dcheck$phase2)
  if (!is.null(postStrata2)){
    ##phase 2 calibration
    for (psvar in postStrata2){
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage==0){
          ## G-calibration at population level
          x2<-qr.resid(psvar$qr,x/psvar$w)*psvar$w
        } else {
          ## G-calibration within clusters
          stop("calibration within clusters not allowed for two-phase designs")
        }
      } else {
        ## ordinary post-stratification
        psw<-attr(psvar, "weights")
        postStrata<-as.factor(psvar)
        psmeans<-rowsum(x/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
        x2<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
      }
    }
  } else x2<-x
  phase2var<-htvar.matrix(x2,design$dcheck$phase2)
  rval<-phase1var+phase2var
  attr(rval,"phases")<-list(phase1=phase1var,phase2=phase2var)
  rval
}

svytotal.twophase2<-function(x,design, na.rm=FALSE, deff=FALSE,...){
    
    
    if (inherits(x,"formula")){
        ## do the right thing with factors
        mf<-model.frame(x,design$phase1$sample$variables,
                        na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
        cols<-sapply(xx,NCOL)
        x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
        scols<-c(0,cumsum(cols))
        for(i in 1:length(xx)){
            x[,scols[i]+1:cols[i]]<-xx[[i]]
        }
        colnames(x)<-do.call("c",lapply(xx,colnames))
    } else {
        if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        else {
            if(is.data.frame(x) && any(sapply(x,is.factor))){
                xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
                cols<-sapply(xx,NCOL)
                scols<-c(0,cumsum(cols))
                cn<-character(sum(cols))
                for(i in 1:length(xx))
                    cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
                x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
                for(i in 1:length(xx)){
                    x[,scols[i]+1:cols[i]]<-xx[[i]]
                }
                colnames(x)<-cn
            }
        }
    }
  x<-as.matrix(x)
  
  if (na.rm){
      nas<-rowSums(is.na(x))
      design<-design[nas==0,]
      if(length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
      else
          x[nas>0,]<-0
  }
  
  N<-sum(1/design$prob)
  total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
  class(total)<-"svystat"
  attr(total, "var")<-v<-twophase2var(x/design$prob,design)
  attr(total,"statistic")<-"total"
  
  if (is.character(deff) || deff){
    nobs<-NROW(design$cluster)
    if (deff=="replace")
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/N
    else
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2
    attr(total, "deff")<-v/vsrs
  }
  
  
  return(total)
}

svymean.twophase2<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$phase1$sample$variables
                    ,na.action=na.pass)
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
  else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else {
          if(is.data.frame(x) && any(sapply(x,is.factor))){
              xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
              cols<-sapply(xx,NCOL)
              scols<-c(0,cumsum(cols))
              cn<-character(sum(cols))
              for(i in 1:length(xx))
                  cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
              x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
              for(i in 1:length(xx)){
                  x[,scols[i]+1:cols[i]]<-xx[[i]]
              }
              colnames(x)<-cn
          }
      }
  }
  
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    if (any(nas>0))
      design<-design[nas==0,]
    x[nas>0,]<-0
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-twophase2var(x*pweights/psum,design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
      nobs<-nrow(design)
      if(deff=="replace"){
        vsrs<-svyvar(x,design,na.rm=na.rm)/(nobs)
      } else {
        if(psum<nobs) {
          vsrs<-NA*v
          warning("Sample size greater than population size: are weights correctly scaled?")
        } else{
          vsrs<-svyvar(x,design,na.rm=na.rm)*(psum-nobs)/(psum*nobs)
        }
      }
      attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}

model.frame.twophase2<-function(formula,phase=2,...){
  if (phase==1)
    formula$phase1$full$variables
  else 
    formula$phase1$sample$variables
}

svyratio.twophase2<-function(numerator=formula, denominator, design, separate=FALSE,na.rm=FALSE,formula,...){

    if (separate){
      strats<-sort(unique(design$phase2$strata[,1]))
      if (!design$phase2$has.strata)
        warning("Separate and combined ratio estimators are the same for unstratified designs")
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$phase2$strata[,1] %in% s),
                                   separate=FALSE,...)
                     attr(tmp,"call")<-bquote(Stratum==.(s))
                     tmp}))
      names(rval$ratios)<-strats
   
      class(rval)<-c("svyratio_separate")
      rval$call<-sys.call()
      rval$strata<-strats
      return(rval)
    }
  
    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,model.frame(design),na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,model.frame(design),na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, model.frame(design))

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    nas<-!complete.cases(all)
    if (na.rm){
        design<-design[!nas,]
        all<-all[!nas,,drop=FALSE]
        numerator<-numerator[!nas,,drop=FALSE]
        denominator<-denominator[!nas,,drop=FALSE]
    }
    allstats<-svytotal(all, design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-twophase2var(r*1/design$prob, design)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }


"[.twophase2"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
    ## Set weights to zero:  don't try to save memory
    ## There should be an easier way to complement a subscript..
    if (is.logical(i) && any(!i)){
      ## logical indexing: use !
      x$prob[!i]<-Inf
      x$phase2$prob[!i]<-Inf
      x$dcheck<-lapply(x$dcheck, function(m) {m[!i,!i]<-0; m})
    } else if (is.numeric(i) && length(i)){
      ## numeric indexing: use -
      x$prob[-i]<-Inf
      x$phase2$prob[-i]<-Inf
      x$dcheck<-lapply(x$dcheck, function(m) {m[-i,-i]<-0;m})
    } else if (is.character(i)){
      ##character indexing: use brute force and ignorance
      tmp<-x$prob[i,]
      x$prob<-rep(Inf, length(x$prob))
      x$prob[i,]<-tmp
      tmp<-x$phase2$prob[i,]
      x$phase2$prob<-rep(Inf, length(x$phase2$prob))
      x$phase2$prob[i,]<-tmp
      x$dcheck<-lapply(x$dcheck, function(m) {n<-Matrix(ncol(m),ncol(m)); n[i,i]<-m[i,i]})
    }
    index<-is.finite(x$prob)
    psu<-!duplicated(x$phase2$cluster[index,1])
    tt<-table(x$phase2$strata[index,1][psu])
    if(any(tt==1)){
      warning(sum(tt==1)," strata have only one PSU in this subset.")
    }
  } else {
    x$phase1$full<-x$phase1$full[,...]
    x$phase1$sample<-x$phase1$sample[,...]
    x$phase2<-x$phase2[,...]
  }
  x
}

dim.twophase2<-function(x,...){
	dim(x$phase1$sample$variables)
}

degf.twophase2<-function(design,...) degf(design$phase2)

na.fail.twophase2<-function(object,...){
	tmp<-na.fail(object$phase1$sample$variables,...)
	object
}

na.omit.twophase2<-function(object,...){
  tmp<-na.omit(object$phase1$sample$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object<-object[-omit,]
    object$phase1$sample$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

na.exclude.twophase2<-function(object,...){
	tmp<-na.exclude(object$phase1$sample$variables,...)
	exclude<-attr(tmp,"na.action")
	if (length(exclude)){
           object<-object[-exclude,]
	   object$phase1$sample$variables<-tmp
	   attr(object,"na.action")<-exclude
	}
	object
}


update.twophase2<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$phase1$sample$variables[,newnames[j]]<-eval(dots[[j]], object$phase1$sample$variables, parent.frame())
    object$phase1$full$variables[,newnames[j]]<-eval(dots[[j]], object$phase1$full$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}

subset.twophase2<-function(x,subset,...){
  e <- substitute(subset)
  r <- eval(e, x$phase1$sample$variables, parent.frame())
  r <- r & !is.na(r) 
  x<-x[r,]
  x$call<-sys.call(-1)
  x
}


calibrate.twophase2<-function(design, phase=2, formula, population,
                             calfun=c("linear","raking","logit","rrz"),...){

    if (phase==1){
        phase1<-calibrate(design$phase1$full,formula, population, ...)
        design$phase1$full<-phase1
        design$phase1$sample<-phase1[design$subset,]
        
    } else if(phase==2){

        if (is.character(calfun)) calfun<-match.arg(calfun)
        if (is.character(calfun) && calfun=="rrz"){
            design<-estWeights(design, formula,...)
            design$call<-sys.call(-1)
            return(design)
        }
            
        if (missing(population) || is.null(population)){
            ## calibrate to phase 1 totals
            population<-colSums(model.matrix(formula,
                               model.frame(formula, design$phase1$full$variables)))
        }
        
        phase2<-design$phase2
        phase2$variables<-design$phase1$sample$variables
        phase2<-calibrate(phase2,formula,population,calfun=calfun,...)
        g<-design$phase2$prob/phase2$prob
        phase2$variables<-NULL
        design$phase2<-phase2
	design$usu<-design$usu/g

    } else stop("`phase' must be 1 or 2")

    
    if (length(design$phase1$sample$prob)==length(design$phase2$prob))
        design$prob<-design$phase1$sample$prob*design$phase2$prob
    else{
        design$prob<-rep(Inf,length(design$phase1$sample$prob))
        design$prob[subset]<-design$prob[subset]*design$phase2$prob
    }

    design$call<-sys.call(-1)
    
    design

}


postStratify.twophase2<-function(design, ...) {
	stop("postStratify not yet implemented for two-phase designs. Use calibrate()")
}

                             
estWeights.twophase2<-function(data, formula=NULL, working.model=NULL,...){

  if (!xor(is.null(formula), is.null(working.model)))
    stop("Must specify one of formula, working.model")

  certainty<-rep(FALSE,nrow(data$phase1$full$variables))
  certainty[data$subset]<-data$phase2$fpc$popsize==data$phase2$fpc$sampsize

  if (!is.null(formula)){
    ff<-data$subset~rhs
    ff[[3]]<-formula[[2]]
    if(!attr(terms(ff),"intercept")) stop("formula must have an intercept")
    
    model<-glm(ff, data=data$phase1$full$variables, family=binomial(),
               subset=!certainty, na.action=na.fail)
  } else {
    xx<-estfun(working.model)
    model<-glm(data$subset~xx,family=binomial(), subset=!certainty, na.action=na.fail)
  }
  fitp<-as.numeric(certainty[data$subset])
  fitp[!certainty[data$subset]]<-fitted(model)[data$subset[!certainty]]
  
  g<- (1/fitp)/(1/data$phase2$prob)
  
  mm<-model.matrix(model)[data$subset[!certainty],,drop=FALSE]

  if (any(certainty)){
    mm1<-matrix(0,ncol=ncol(mm)+1,nrow=sum(data$subset))
    mm1[,1]<-as.numeric(certainty[data$subset])
    mm1[!certainty[data$subset],-1]<-mm
        mm<-mm1
  }
  
  whalf<-sqrt(1/data$phase2$prob)
  
  caldata<-list(qr=qr(mm*whalf), w=g*whalf, stage=0, index=NULL)
  class(caldata) <- c("greg_calibration","gen_raking")
  
  data$phase2$prob<-fitp
  data$usu<-data$usu/g
  data$phase2$postStrata <- c(data$phase2$postStrata, list(caldata))
    
  if (length(data$phase1$sample$prob)==length(data$phase2$prob))
    data$prob<-data$phase1$sample$prob*data$phase2$prob
  else{
    data$prob<-rep(Inf,length(data$phase1$sample$prob))
    data$prob[subset]<-data$prob[subset]*data$phase2$prob
  }
  
  data$call <- sys.call(-1)
  
  data
  
}


