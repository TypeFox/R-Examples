##
##
twophase<-function(id,strata=NULL, probs=NULL, weights=NULL, fpc=NULL,
                   subset, data, method=c("full","approx","simple")){
  method<-match.arg(method)
  if(method=="full") {
    if (!is.null(weights)) stop("weights not accepted by method='full'")
    return(twophase2(id=id, strata=strata, probs=probs, fpc=fpc,subset=subset,data=data))
  }
                     
  d1<-svydesign(ids=id[[1]],strata=strata[[1]],weights=weights[[1]],
                probs=probs[[1]],fpc=fpc[[1]],data=data)

  if(inherits(subset,"formula"))
    subset<-eval.parent(model.frame(subset,data=data,na.action=na.pass))[[1]]

  if(!is.logical(subset) && sort(unique(subset))==c(0,1))
      subset<-as.logical(subset)

  if (any(is.na(subset))) stop("missing values in 'subset'")
  
  d1s<-svydesign(ids=id[[1]],strata=strata[[1]],weights=weights[[1]],
                probs=probs[[1]],fpc=fpc[[1]],data=data[subset,])
  d1s$prob<-d1$prob[subset]
  d1s$allprob<-d1$allprob[subset,,drop=FALSE]

  ##if (NCOL(d1s$allprob)>1)
  ##  stop("Can't handle multistage sampling at phase 1 (yet)")
  
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
                weights=weights[[2]], fpc=fpc[[2]], data=data[subset,])

  ## ugly hack to get nicer labels
  if(!is.null(fpc[[2]])){
    d2call<-bquote(svydesign(ids=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                              weights=.(weights[[2]]), fpc=.(fpc[[2]])))
  } else{
    d2call<-bquote(svydesign(ids=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                              weights=.(weights[[2]]), fpc=`*phase1*`))
  }
  for(i in names(d2call)[-1])
    d2call[[i]]<-d2call[[i]]
  d2$call<-d2call
  d1call<-bquote(svydesign(ids=.(id[[1]]), strata=.(strata[[1]]), probs=.(probs[[1]]),
                              weights=.(weights[[1]]), fpc=.(fpc[[1]])))
  for(i in names(d1call)[-1])
    d1call[[i]]<-d1call[[i]]
  d1$call<-d1call


  ## Add phase 2 fpc and probs if they were computed rather than specified.
  if (!is.null(popsize))
    d2$fpc<-as.fpc(popsize[subset,,drop=FALSE],d2$strata,d2$cluster)
  if(is.null(probs[[2]]) && is.null(weights[[2]]) && !is.null(d2$fpc$popsize)){
    d2$allprob<-1/weights(d2$fpc,final=FALSE)
    d2$prob<-apply(as.data.frame(d2$allprob),1,prod)
  }
  
  d2$variables<-NULL
  rval<-list(phase1=list(full=d1,sample=d1s),
             phase2=d2,
             subset=subset)
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

##  if (any(rval$usu<1) && any(duplicated(d1$cluster[,1])))
##      stop("Phase 1 design must either be element sampling or have all phase 1 sampling units in phase 2")
  
  if (length(rval$phase1$sample$prob)==length(d2$prob))
    rval$prob<-rval$phase1$sample$prob*d2$prob
  else{
    rval$prob<-rep(Inf,length(rval$phase1$sample$prob))
    rval$prob[subset]<-rval$prob[subset]*d2$prob
  }
  rval$call<-sys.call()
  class(rval) <- c("twophase","survey.design")
  rval
}

print.twophase<-function(x,...){
  cat("Two-phase design: ")
  print(x$call)
  cat("Phase 1:\n")
  print(x$phase1$full)
  cat("Phase 2:\n")
  print(x$phase2)
  invisible(x)
}

summary.twophase<-function(object,...){
  class(object)<-"summary.twophase"
  object
}

print.summary.twophase<-function(x,...,varnames=TRUE){
  cat("Two-phase design: ")
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

twophasevar<-function(x,design){
  d1 <- design$phase1$sample
  if (NROW(x)==length(design$usu)){
      ph2pr<-design$usu
      if (any(design$prob==Inf))
          x[is.na(x)]<-0
  }else{
      x[is.na(x)]<-0
      ph2pr<-rep(1,NROW(x))
      ph2pr[design$subset]<-design$usu
  }
  ## compute phase 1 variance
  vphase1 <- svyrecvar.phase1(x,d1$cluster, d1$strata, d1$fpc,
                              postStrata=d1$postStrata,
                              ph2prob=ph2pr,
                              nPSUfull=design$phase1$full$fpc$sampsize[design$subset,,drop=FALSE])
  
  ## is phase 2 sampling whole phase 1 units or subsampling within units?
  if (design$samescale)
    u2<-x
  else
    u2<-x*sqrt(d1$prob)

  u2[is.na(u2)]<-0
  ## compute phase 2 variance
  vphase2 <- with(design, svyrecvar(u2, phase2$cluster, phase2$strata,
                                  phase2$fpc, postStrata=phase2$postStrata))
  rval <- vphase1+vphase2
  attr(rval, "phases")<-list(phase1=vphase1, phase2=vphase2)
  rval
}

svyrecvar.phase1<-function(x, clusters,  stratas, fpcs, postStrata=NULL,
                           lonely.psu=getOption("survey.lonely.psu"),
                           one.stage=getOption("survey.ultimate.cluster"),
                           ph2prob, nPSUfull){
    
    x<-as.matrix(x)
    cal<-NULL

    ## FIXME: calibration of phase 1 not yet implemented.
    ## Remove post-stratum means, which may cut across clusters
    ## Also center the data using any "g-calibration" models
    if(!is.null(postStrata)){
        stop("calibration of phase 1 not yet implemented")
        for (psvar in postStrata){
            if (inherits(psvar, "greg_calibration")) {
                if (psvar$stage==0){
                    ## G-calibration at population level
                    x<-qr.resid(psvar$qr,x/psvar$w)*psvar$w
                } else {
                    ## G-calibration within clusters
                    cal<-c(cal, list(psvar))
                }
            } else {
                ## ordinary post-stratification
                psw<-attr(psvar, "weights")
                postStrata<-as.factor(psvar)
                psmeans<-rowsum(x/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
                x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
            }
        }
    }
  
    multistage.phase1(x, clusters,stratas,fpcs$sampsize, fpcs$popsize,
             lonely.psu=getOption("survey.lonely.psu"),
             one.stage=one.stage,stage=1,cal=cal,ph2prob=ph2prob,
                      nPSUfull=nPSUfull)
}


multistage.phase1<-function(x, clusters,  stratas, nPSUs, fpcs,
                    lonely.psu=getOption("survey.lonely.psu"),
                     one.stage=FALSE,stage,cal,ph2prob, nPSUfull){
  
  n<-NROW(x)
 
  
  v <- onestage.phase1(x,stratas[,1], clusters[,1], nPSUs[,1],
                fpcs[,1], lonely.psu=lonely.psu,stage=stage,cal=cal,
                ph2prob=ph2prob, nPSUfull=nPSUfull[,1])
  
  if (one.stage!=TRUE && !is.null(fpcs) && NCOL(clusters)>1) {
    v.sub<-by(1:n, list(as.numeric(clusters[,1])), function(index){
      ## residuals for G-calibration using population information
      ## only on clusters at this stage.
      for(cali in cal){
        if (cali$stage != stage)
          next
        j<-match(clusters[index,1],cali$index)
        if (length(unique(j))!=1)
          stop("Internal problem in g-calibration data: stage",stage,
               ", cluster", j)
        j<-j[[1]]
        x[index,]<-qr.resid(cali$qr[[j]], x[index,,drop=FALSE]/cali$w[[j]])*cali$w[[j]]
      }
      multistage.phase1(x[index,,drop=FALSE], clusters[index,-1,drop=FALSE],
                        stratas[index,-1,drop=FALSE], nPSUs[index,-1,drop=FALSE],
                        fpcs[index,-1,drop=FALSE],
                        lonely.psu=lonely.psu,one.stage=one.stage-1,
                        stage=stage+1,cal=cal,ph2prob=ph2prob[index],
                        nPSUfull=nPSUfull[index,-1,drop=FALSE])*nPSUfull[index[1],1]/fpcs[index[1],1]
    })
    
    for(i in 1:length(v.sub))
      v<-v+v.sub[[i]]
  }
  v
}


onestrat.phase1<-function(x,cluster,nPSU,fpc, lonely.psu,stratum=NULL,
                          stage=1,cal,ph2prob, nPSUfull){
  x<-rowsum(x, cluster)
  ph2prob<-ph2prob[!duplicated(cluster)]

  nsubset<-nrow(x)
  if (nsubset<nPSU)
     x<-rbind(x,matrix(0,ncol=ncol(x),nrow=nPSU-nrow(x)))

  ph2prob<-c(ph2prob,rep(1,nPSU-nsubset))
  xcenter<-colMeans(x*nPSU/nPSUfull)

  x<-x*ph2prob

  
  if (is.null(fpc))
      f<-1
  else
      f<-ifelse(fpc==Inf, 1, (fpc-nPSUfull)/fpc)
  
  if (lonely.psu!="adjust" || nsubset>1 ||
      (nPSU>1 && !getOption("survey.adjust.domain.lonely")))
      x<-sweep(x, 2, xcenter, "-")
  
  if (nPSU>1)
      scale<-f*nPSUfull/(nPSUfull-1)
  else
      scale<-f

  if (nsubset==1 && nPSU>1){ 
      warning("Stratum (",stratum,") has only one PSU at stage ",stage)
      if (lonely.psu=="average" && getOption("survey.adjust.domain.lonely"))
          scale<-NA
    }
  
  if (nPSU>1){
      return(crossprod(x/sqrt(ph2prob))*scale)
  } else if (f<0.0000001) ## certainty PSU
      return(0*crossprod(x/sqrt(ph2prob)))
  else {
      rval<-switch(lonely.psu,
                   certainty=scale*crossprod(x/sqrt(ph2prob)),
                   remove=scale*crossprod(x/sqrt(ph2prob)),
                   adjust=scale*crossprod(x/sqrt(ph2prob)),
                   average=NA*crossprod(x/sqrt(ph2prob)),
                   fail= stop("Stratum (",stratum,") has only one PSU at stage ",stage),
                   stop("Can't handle lonely.psu=",lonely.psu)
            )
      rval
  }
}


onestage.phase1<-function(x, strata, clusters, nPSU, fpc,
                          lonely.psu=getOption("survey.lonely.psu"),stage=0,
                          cal,ph2prob, nPSUfull){
  stratvars<-tapply(1:NROW(x), list(factor(strata)), function(index){
    onestrat.phase1(x[index,,drop=FALSE], clusters[index],
                    nPSU[index][1], fpc[index][1],
                    lonely.psu=lonely.psu,stratum=strata[index][1], stage=stage,cal=cal,
                    ph2prob=ph2prob[index], nPSUfull=nPSUfull[index][1])
  })
  p<-NCOL(x)
  nstrat<-length(unique(strata))
  nokstrat<-sum(sapply(stratvars,function(m) !any(is.na(m))))
  apply(array(unlist(stratvars),c(p,p,length(stratvars))),1:2,sum,na.rm=TRUE)*nstrat/nokstrat
}


svytotal.twophase<-function(x,design, na.rm=FALSE, deff=FALSE,...){
    
    
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
  attr(total, "var")<-v<-twophasevar(x/design$prob,design)
  attr(total,"statistic")<-"total"
  
  if (is.character(deff) || deff){
    nobs<-NROW(design$cluster)
    if (deff=="replace")
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design)^2)*(N-nobs)/N
    else
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design)^2)
    attr(total, "deff")<-v/vsrs
  }
  
  
  return(total)
}

svymean.twophase<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
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
    design<-design[nas==0,]
    if(length(nas)>length(design$prob))
        x<-x[nas==0,,drop=FALSE]
    else
        x[nas>0,]<-0
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-twophasevar(x*pweights/psum,design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
      nobs<-NROW(design$cluster)
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

model.frame.twophase<-function(formula,phase=2,...){
  if (phase==1)
    formula$phase1$full$variables
  else 
    formula$phase1$sample$variables
}

svyratio.twophase<-function(numerator=formula, denominator, design, separate=FALSE,na.rm=FALSE,formula,...){

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
        vars[i,j]<-twophasevar(r*1/design$prob, design)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }


"[.twophase"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
      if (is.calibrated(x$phase1$full) || is.calibrated(x$phase2) || !drop){
          ## Set weights to zero: no memory saving possible
          ## There should be an easier way to complement a subscript..
          if (is.logical(i)){
              x$prob[!i]<-Inf
              x$phase2$prob[!i]<-Inf
          } else if (is.numeric(i) && length(i)){
              x$prob[-i]<-Inf
              x$phase2$prob[-i]<-Inf
          } else {
              tmp<-x$prob[i,]
              x$prob<-rep(Inf, length(x$prob))
              x$prob[i,]<-tmp
          }
          index<-is.finite(x$prob)
          psu<-!duplicated(x$phase2$cluster[index,1])
          tt<-table(x$phase2$strata[index,1][psu])
          if(any(tt==1)){
              warning(sum(tt==1)," strata have only one PSU in this subset.")
          }
      } else {
          ## subset everything.
          x$prob<-x$prob[i]
          if (is.logical(i))
              x$subset[x$subset]<- i
          else if (is.numeric(i) && length(i))
              x$subset[which(x$subset)[-i]]<- FALSE
          else
              x$subset<-FALSE & x$subset
          x$usu<-x$usu[i]
          x$phase1$sample<-x$phase1$sample[i,...,drop=TRUE]
          x$phase2<-x$phase2[i,...,drop=TRUE]
        }
  } else {
    x$phase1$full<-x$phase1$full[,...]
    x$phase1$sample<-x$phase1$sample[,...]
    x$phase2<-x$phase2[,...]
  }
  x
}

dim.twophase<-function(x,...){
	dim(x$phase1$sample$variables)
}

na.fail.twophase<-function(object,...){
	tmp<-na.fail(object$phase1$sample$variables,...)
	object
}

na.omit.twophase<-function(object,...){
  tmp<-na.omit(object$phase1$sample$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object<-object[-omit,]
    object$phase1$sample$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

na.exclude.twophase<-function(object,...){
	tmp<-na.exclude(object$phase1$sample$variables,...)
	exclude<-attr(tmp,"na.action")
	if (length(exclude)){
           object<-object[-exclude,]
	   object$phase1$sample$variables<-tmp
	   attr(object,"na.action")<-exclude
	}
	object
}


update.twophase<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$phase1$sample$variables[,newnames[j]]<-eval(dots[[j]], object$phase1$sample$variables, parent.frame())
    object$phase1$full$variables[,newnames[j]]<-eval(dots[[j]], object$phase1$full$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}

subset.twophase<-function(x,subset,...){
        e <- substitute(subset)
        r <- eval(e, x$phase1$sample$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
	x$call<-sys.call(-1)
	x
}


calibrate.twophase<-function(design, phase=2, formula, population,
                             calfun=c("linear","raking","logit","rrz"),...){

    if (phase==1){
        stop("phase 1 calibration not yet implemented")
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


postStratify.twophase<-function(design, ...) {
	stop("postStratify not implemented for two-phase designs. Use calibrate()")
}

estWeights<-function(data, formula, ...) UseMethod("estWeights")
                             
estWeights.twophase<-function(data, formula=NULL, working.model=NULL,...){

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


estfun<-function(model,...) UseMethod("estfun")
estfun.coxph<-function(model, ...) resid(model,"score")
estfun.glm<-function(model){
  xmat<-model.matrix(model)
  residuals(model,"working")*model$weights*xmat
}
estfun.lm<-function(model,...){
  model.matrix(model)*resid(model)
}




estWeights.data.frame<-function(data,formula=NULL, working.model=NULL,
                                subset=NULL, strata=NULL,...){

    if (is.null(subset)){
        subset<-complete.cases(data)
        if (all(subset))
          stop("No missing data.")
        }
    
    if(is.null(strata)){
        des<-twophase(id=list(~1,~1), subset=subset, data=data)
    } else{
        des<-twophase(id=list(~1,~1), subset=subset, data=data,
                      strata=list(NULL,strata))
    }

    rval<-estWeights(des,formula=formula,working.model=working.model)
    rval$call<-sys.call(-1)
    rval
    
}
