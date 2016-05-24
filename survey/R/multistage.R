##
##  Recursive estimation of linearisation variances
##  in multistage samples.
##

svydesign<-function(ids, probs = NULL, strata = NULL, variables = NULL, 
    fpc = NULL, data=NULL, nest = FALSE, check.strata = !nest, 
    weights = NULL,pps=FALSE,...){
	UseMethod("svydesign", data)
	}

svydesign.default<-function(ids,probs=NULL,strata=NULL,variables=NULL, fpc=NULL,
                    data=NULL, nest=FALSE, check.strata=!nest,weights=NULL,pps=FALSE,
                            variance=c("HT","YG"),...){
  variance<-match.arg(variance)
  if(is.character(pps)){
    a<-match.arg(pps,c("brewer","overton","other"))
    if (!(pps %in% c("brewer","other")))
      return(pps_design(ids=ids,probs=probs, strata=strata,variables=variables, fpc=fpc,
                 data=data,method=a,call=sys.call(-1),variance=variance,...))
  } else if (!is.logical(pps)){
    return(pps_design(ids=ids,probs=probs, strata=strata,variables=variables, fpc=fpc,
                      data=data,method=pps,call=sys.call(-1),variance=variance,...))
  }

  if (!is.character(pps) || pps!="other"){
    if (variance!="HT")
      stop("Only variance='HT' supported for this design")
  }
  
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

    na.failsafe<-function(message="missing values in object"){
      function(object,...){
        if (NCOL(object)==0)
          object
        else {
          ok <- complete.cases(object)
          if (all(ok)) 
            object
          else stop(message)
        }
      }
    }

     na.id<-na.failsafe("missing values in `id'")
     if(inherits(ids,"formula")) {
	 mf<-substitute(model.frame(ids,data=data, na.action=na.id))
	 ids<-eval.parent(mf)
         if (ncol(ids)==0) ## formula was ~1
           ids<-data.frame(id=1:nrow(ids))
       } else{
         if (is.null(ids))
           stop("Must provide ids= argument")
         else
           ids<-na.id(data.frame(ids))
       }

    na.prob<-na.failsafe("missing values in `prob'")
    if(inherits(probs,"formula")){
      mf<-substitute(model.frame(probs,data=data,na.action=na.prob))
      probs<-eval.parent(mf)
    }

    na.weight<-na.failsafe("missing values in `weights'")
    if(inherits(weights,"formula")){
      mf<-substitute(model.frame(weights,data=data,na.action=na.weight))
      weights<-eval.parent(mf)
     } else if (!is.null(weights))
         weights<-na.weight(data.frame(weights))
    if(!is.null(weights)){
      if (!is.null(probs))
         stop("Can't specify both sampling weights and probabilities")
       else
         probs<-as.data.frame(1/as.matrix(weights))
     }

      

    na.strata<-na.failsafe("missing values in `strata'")
    if (!is.null(strata)){
      if(inherits(strata,"formula")){
        mf<-substitute(model.frame(strata,data=data, na.action=na.strata))
        strata<-eval.parent(mf)
      }
      if (!is.list(strata))
        strata<-data.frame(strata=strata)
      has.strata<-TRUE
    } else {
      has.strata <-FALSE
      strata<-na.strata(as.data.frame(matrix(1, nrow=NROW(ids), ncol=NCOL(ids))))
    }

    
    if (inherits(variables,"formula")){
        mf<-substitute(model.frame(variables,data=data,na.action=na.pass))
        variables <- eval.parent(mf)
    } else if (is.null(variables)){
        variables<-data
    } else
        variables<-do.call("data.frame",variables)


    na.fpc<-na.failsafe("missing values in `fpc'")
    if (inherits(fpc,"formula")){
      mf<-substitute(model.frame(fpc,data=data,na.action=na.fpc))
      fpc<-eval.parent(mf)
    }
      
  ## check for only one PSU: probably a typo
  if ((length(unique(ids[,1]))==1) && !(nest && has.strata)){
    stop("Design has only one primary sampling unit")
  }
  
      ## force subclusters nested in clusters
      if (NCOL(ids)>1){
        N<-ncol(ids)
        for(i in 2:N){
          ids[,i]<-do.call("interaction", ids[,1:i,drop=FALSE])
        }
      }
      ## force clusters nested in strata
      if (nest && has.strata && NCOL(ids)){
        N<-NCOL(ids)
        NS<-NCOL(strata)
        for(i in 1:N)
          ids[,i]<-do.call("interaction",
                           c(strata[,1:min(i,NS),drop=FALSE], ids[,i,drop=FALSE]))
      }
      
    ## check if clusters nested in strata 
     if (check.strata && nest)
      warning("No point in check.strata=TRUE if nest=TRUE")
    if(check.strata && !is.null(strata) && NCOL(ids)){
       sc<-(rowSums(table(ids[,1],strata[,1])>0))
       if(any(sc>1)) stop("Clusters not nested in strata at top level; you may want nest=TRUE.")
    }

      ## force substrata nested in clusters
      N<-ncol(ids)
      NS<-ncol(strata)
      if (N>1){
        for(i in 2:N)
          strata[,i]<-interaction(strata[,min(i,NS)], ids[,i-1])
      }

    ## PPS: valid choices currently are FALSE and "brewer"
    if (is.logical(pps) && pps) stop("'pps' must be FALSE or a character string")
    if (is.character(pps)) {
      pps<-TRUE
    }
    
    ## Finite population correction: specified per observation
    ## Also incorporates design sample sizes formerly in nPSU
    
      if (!is.null(fpc) && !is.numeric(fpc) && !is.data.frame(fpc))
        stop("fpc must be a matrix or dataframe or NULL")

      fpc<-as.fpc(fpc,strata, ids, pps=pps)

      ## if FPC specified, but no weights, use it for weights
    if (is.null(probs) && is.null(weights)){
      if (is.null(fpc$popsize)){
        if (missing(probs) && missing(weights))
          warning("No weights or probabilities supplied, assuming equal probability")
        probs<-rep(1,nrow(ids))
      } else {
        probs<-1/weights(fpc, final=FALSE)
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
    rval$call<-sys.call(-1)
    rval$pps<-pps
    class(rval)<-c("survey.design2","survey.design")
    rval
}

onestrat<-function(x,cluster,nPSU,fpc, lonely.psu,stratum=NULL,stage=1,cal=cal){
  
  if (is.null(fpc))
      f<-rep(1,NROW(x))
  else{
      f<-ifelse(fpc==Inf, 1, (fpc-nPSU)/fpc)
  }

  if (nPSU>1)
      scale<-f*nPSU/(nPSU-1)
  else
      scale<-f
  if (all(f<0.0000001))## self-representing stratum
      return(matrix(0,NCOL(x),NCOL(x)))

  scale<-scale[!duplicated(cluster)]
  
  x<-rowsum(x,cluster)
  nsubset<-nrow(x)
  
  if (nsubset<nPSU) {
    ##can't be PPS, so scale must be a constant
    x<-rbind(x,matrix(0,ncol=ncol(x),nrow=nPSU-nrow(x)))
    scale<-rep(scale[1],NROW(x))
  }
  if (lonely.psu!="adjust" || nsubset>1 ||
      (nPSU>1 & !getOption("survey.adjust.domain.lonely")))
      x<-sweep(x, 2, colMeans(x), "-")

  if (nsubset==1 && nPSU>1 && getOption("survey.adjust.domain.lonely")){ 
      warning("Stratum (",stratum,") has only one PSU at stage ",stage)
      if (lonely.psu=="average" && getOption("survey.adjust.domain.lonely"))
          scale<-NA
    }
  
  if (nPSU>1){
      return(crossprod(x*sqrt(scale)))
  } else {
      rval<-switch(lonely.psu, 
                   certainty=crossprod(x*sqrt(scale)),
                   remove=crossprod(x*sqrt(scale)),
                   adjust=crossprod(x*sqrt(scale)),
                   average=NA*crossprod(x),
                   fail= stop("Stratum (",stratum,") has only one PSU at stage ",stage),
                   stop("Can't handle lonely.psu=",lonely.psu)
            )
      rval
  }
}


onestage<-function(x, strata, clusters, nPSU, fpc, lonely.psu=getOption("survey.lonely.psu"),stage=0, cal){
  stratvars<- tapply(1:NROW(x), list(factor(strata)), function(index){
             onestrat(x[index,,drop=FALSE], clusters[index],
             nPSU[index][1], fpc[index], ##changed from fpc[index][1], to allow pps(brewer)
             lonely.psu=lonely.psu,stratum=strata[index][1], stage=stage,cal=cal)
  })
  p<-NCOL(x)
  nstrat<-length(unique(strata))
  nokstrat<-sum(sapply(stratvars,function(m) !any(is.na(m))))
  apply(array(unlist(stratvars),c(p,p,length(stratvars))),1:2,sum,na.rm=TRUE)*nstrat/nokstrat
}


svyrecvar<-function(x, clusters,  stratas, fpcs, postStrata=NULL,
                    lonely.psu=getOption("survey.lonely.psu"),
                    one.stage=getOption("survey.ultimate.cluster")){

  x<-as.matrix(x)
  cal<-NULL
  
  ## Remove post-stratum means, which may cut across clusters
  ## Also center the data using any "g-calibration" models
  if(!is.null(postStrata)){
    for (psvar in postStrata){
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage==0){
          ## G-calibration at population level
          x<-qr.resid(psvar$qr,x/psvar$w)*psvar$w
        } else {
          ## G-calibration within clusters
          cal<-c(cal, list(psvar))
        }
      } else if (inherits(psvar, "raking")){
        ## raking by iterative proportional fitting
        for(iterations in 1:10){
          for(margin in psvar){
            psw<-attr(margin, "weights")
            x<- x - psw*apply(x/psw, 2, ave, margin)
          }
        }
      } else {
        ## ordinary post-stratification
        psw<-attr(psvar, "weights")
        oldw<-attr(psvar, "oldweights")
        if (is.null(oldw)) oldw<-rep(1,length(psw))
        zeroes<-which(psw==0 & oldw==0)
        if (length(zeroes)) psw[zeroes]=1
        psvar<-as.factor(psvar)
        psmeans<-rowsum(x*oldw/psw,psvar,reorder=TRUE)/as.vector(by(oldw,psvar,sum))
        x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
      }
    }
  }
  
  multistage(x, clusters,stratas,fpcs$sampsize, fpcs$popsize,
             lonely.psu=getOption("survey.lonely.psu"),
             one.stage=one.stage,stage=1,cal=cal)
}

multistage<-function(x, clusters,  stratas, nPSUs, fpcs,
                    lonely.psu=getOption("survey.lonely.psu"),
                     one.stage=FALSE,stage,cal){
  
  n<-NROW(x)
 
  
  v <- onestage(x,stratas[,1], clusters[,1], nPSUs[,1],
                fpcs[,1], lonely.psu=lonely.psu,stage=stage,cal=cal)
  
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
      multistage(x[index,,drop=FALSE], clusters[index,-1,drop=FALSE],
                 stratas[index,-1,drop=FALSE], nPSUs[index,-1,drop=FALSE],
                 fpcs[index,-1,drop=FALSE],
                 lonely.psu=lonely.psu,one.stage=one.stage-1,
                 stage=stage+1,cal=cal)*nPSUs[index[1],1]/fpcs[index[1],1]
    })
    
    for(i in 1:length(v.sub))
      v<-v+v.sub[[i]]
  }
  dimnames(v)<-list(colnames(x),colnames(x))
  v
}


## fpc not given are zero: full sampling.
as.fpc<-function(df,strata,ids,pps=FALSE){

  count<-function(x) length(unique(x))
  
  sampsize<-matrix(ncol=ncol(ids),nrow=nrow(ids))
  for(i in 1:ncol(ids))
    split(sampsize[,i],strata[,i])<-lapply(split(ids[,i],strata[,i]),count)
  
  if (is.null(df)){
    ## No fpc
    rval<-list(popsize=NULL, sampsize=sampsize)
    class(rval)="survey_fpc"
    return(rval)
  }
  
  fpc<-as.matrix(df)
  if (xor(ispopsize<-any(df>1), all(df>=1))){
    big<-which(fpc>=1,arr.ind=TRUE)
    small<-which(fpc<1,arr.ind=TRUE)
    cat("record",big[1,1]," stage",big[1,2],": fpc=", fpc[big[1,,drop=FALSE]],"\n")
    cat("record",small[1,1]," stage ",small[1,2],": fpc=", fpc[small[1,,drop=FALSE]],"\n")      
    stop("Must have all fpc>=1 or all fpc<=1")
  }
  
  if (ispopsize){
    if(pps) stop("fpc must be specified as sampling fraction for PPS sampling")
    popsize<-fpc
  } else {
    popsize<-sampsize/(fpc)
  }
  if (any(popsize<sampsize)){
    toobig<-which(popsize<sampsize,arr.ind=TRUE)
    cat("record",toobig[1,1],"stage",toobig[1,2],": popsize=",popsize[toobig[1,,drop=FALSE]],
        " sampsize=", sampsize[toobig[1,,drop=FALSE]],"\n")
    stop("FPC implies >100% sampling in some strata")
  }
  if (!ispopsize && any(is.finite(popsize) & (popsize>1e10))){
    big<-which(popsize>1e10 & is.finite(popsize),arr.ind=TRUE)
    warning("FPC implies population larger than ten billion (record",big[1,1]," stage ",big[1,2],")")
  }
  if(!pps){
    ## check that fpc is constant within strata.
    for(i in 1:ncol(popsize)){
      diff<-by(popsize[,i], list(strata[,i]), count)
      if (any(as.vector(diff)>1)){
        j<-which(as.vector(diff)>1)[1]
        warning("`fpc' varies within strata: stratum ",names(diff)[j], " at stage ",i)
      }
    }
  } else{
    ## check that fpc is constant with clusters
     diff<-by(popsize[,i], list(ids[,i]), count)
      if (any(as.vector(diff)>1)){
        j<-which(as.vector(diff)>1)[1]
        warning("`fpc' varies within cluster: cluster ",names(diff)[j], " at stage ",i)
      }
   }
  
  
  rval<-list(popsize=popsize, sampsize=sampsize)
  class(rval)<-"survey_fpc"
  rval
}

"weights.survey_fpc"<-function(object,final=TRUE,...){
  if (is.null(object$popsize) || any(object$popsize>1e12))
    stop("Weights not supplied and can't be computed from fpc.")
  if (final) {
    pop<-apply(object$popsize,1,prod)
    samp<-apply(object$sampsize,1,prod)
    pop/samp
  } else {
    object$popsize/object$sampsize
  }
}


    

print.survey.design2<-function(x,varnames=FALSE,design.summaries=FALSE,...){
  n<-NROW(x$cluster)
  if (x$has.strata) cat("Stratified ")
  un<-length(unique(x$cluster[,1]))
  if(n==un){
    cat("Independent Sampling design")
    is.independent<-TRUE
    if (is.null(x$fpc$popsize))
      cat(" (with replacement)\n")
    else cat("\n")
  } else {
    cat(NCOL(x$cluster),"- level Cluster Sampling design")
    if (is.null(x$fpc$popsize))
      cat(" (with replacement)\n")
    else cat("\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    cat(paste("With (",paste(unlist(nn),collapse=", "),") clusters.\n",sep=""))
    is.independent<-FALSE
  }

  print(x$call)
  if (design.summaries){
    cat("Probabilities:\n")
    print(summary(x$prob))
    if(x$has.strata){
      if (NCOL(x$cluster)>1)
        cat("First-level ")
      cat("Stratum Sizes: \n")
      oo<-order(unique(x$strata[,1]))
      a<-rbind(obs=table(x$strata[,1]),
	       design.PSU=x$fpc$sampsize[!duplicated(x$strata[,1]),1][oo],
               actual.PSU=table(x$strata[!duplicated(x$cluster[,1]),1]))
      print(a)
    }
    if (!is.null(x$fpc$popsize)){
      if (x$has.strata) {
        cat("Population stratum sizes (PSUs): \n")
        s<-!duplicated(x$strata[,1])
        a<-x$fpc$popsize[s,1]
        names(a)<-x$strata[s,1]
        a<-a[order(names(a))]
        print(a)
      } else {
        cat("Population size (PSUs):",x$fpc$popsize[1,1],"\n")
      }
    }
  }
  if (varnames){
    cat("Data variables:\n")
    print(colnames(x))
  }
  invisible(x)
}

    
summary.survey.design2<-function(object,...){
  class(object)<-c("summary.survey.design2",class(object))
  object
}

print.summary.survey.design2<-function(x,...){
  y<-x
  class(y)<-c("survey.design2",class(x))
  print(y,varnames=TRUE,design.summaries=TRUE,...)
}	
     

.svycheck<-function(object){
  if (inherits(object,"survey.design") &&
      !is.null(object$nPSU))
    warning("This is an old-style design object. Please use as.svydesign2 to update it.")
}

as.svydesign2<-function(object){
  if (inherits(object,"survey.design2"))
    return(object)
  if (!inherits(object,"survey.design"))
    stop("This function is for updating old-style survey.design objects")
  

  count<-function(x) length(unique(x))
  
  strata<-data.frame(one=object$strata)
  if ((nc<-ncol(object$cluster))>1){
    for(i in 2:nc){
      strata<-cbind(strata,object$cluster[,i-1])
    }
  }
  
  sampsize<-matrix(ncol=nc,nrow=nrow(object$cluster))
  
  sampsize[,1]<-object$nPSU[match(object$strata, names(object$nPSU))]
  if (nc>1){
    for(i in 2:nc){
      split(sampsize[,i],strata[,i])<-lapply(split(object$cluster[,i],strata[,i]),count)
    }
  }
  
  if (!is.null(object$fpc)){
    popsize<-sampsize
    popsize[,1]<-object$fpc$N[match(object$strata,object$fpc$strata)]
  } else popsize<-NULL
  if (nc>1 && !is.null(object$fpc)){
    warning("Assuming complete sampling at stages 2 -",nc)
  }

  fpc<-list(popsize=popsize,sampsize=sampsize)
  class(fpc)<-"survey_fpc"
  
           
  object$fpc<-fpc
  object$strata<-strata
  object$nPSU<-NULL
  class(object)<-c("survey.design2","survey.design")
  object
  
}

is.pps<-function(x) if(is.null(x$pps)) FALSE else (x$pps!=FALSE)
    
"[.survey.design2"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
      if (is.calibrated(x) || is.pps(x) || !drop){
          ## Set weights to zero: no memory saving possible
          ## There should be an easier way to complement a subscript..
          if (is.logical(i))
              x$prob[!i]<-Inf
          else if (is.numeric(i) && length(i))
              x$prob[-i]<-Inf
          else {
              tmp<-x$prob[i,]
              x$prob<-rep(Inf, length(x$prob))
              x$prob[i,]<-tmp
          }
          index<-is.finite(x$prob)
          psu<-!duplicated(x$cluster[index,1])
          tt<-table(x$strata[index,1][psu])
          if(any(tt==1) && getOption("survey.adjust.domain.lonely")){
              warning(sum(tt==1)," strata have only one PSU in this subset.")
          }
      } else {
          ## subset everything.
          if (!is.null(x$variables)) ## phase 2 of twophase design
              x$variables<-"[.data.frame"(x$variables,i,..1,drop=FALSE)
          x$cluster<-x$cluster[i,,drop=FALSE]
          x$prob<-x$prob[i]
          x$allprob<-x$allprob[i,,drop=FALSE]
          x$strata<-x$strata[i,,drop=FALSE]
          x$fpc$sampsize<-x$fpc$sampsize[i,,drop=FALSE]
          x$fpc$popsize<-x$fpc$popsize[i,,drop=FALSE]
      }
      
  } else {
      if(!is.null(x$variables))
          x$variables<-x$variables[,..1,drop=FALSE]
  }
  
  x
}

svytotal.survey.design2<-function(x,design, na.rm=FALSE, deff=FALSE,...){

  
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
    } else{
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
        if (length(nas)>length(design$prob))
            x<-x[nas==0,,drop=FALSE]
        else
            x[nas>0,]<-0
    }

    N<-sum(1/design$prob)
    total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
    class(total)<-"svystat"
    attr(total, "var")<-v<-svyrecvar(x/design$prob,design$cluster,
                                     design$strata, design$fpc,
                                   postStrata=design$postStrata)
    attr(total,"statistic")<-"total"

    if (is.character(deff) || deff){
      nobs<-sum(weights(design)!=0)
      if (deff=="replace")
        vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2/nobs
      else
        vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/(N*nobs)
      attr(total, "deff")<-v/vsrs
    }
    

  return(total)
}


svymean.survey.design2<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
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
  else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else if(is.data.frame(x) && any(sapply(x,is.factor))){
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
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    if (length(nas)>length(design$prob))
        x<-x[nas==0,,drop=FALSE]
    else
        x[nas>0,]<-0
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-svyrecvar(x*pweights/psum,design$cluster,design$strata, design$fpc,
              postStrata=design$postStrata)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
      nobs<-sum(weights(design)!=0)
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

svyratio.survey.design2<-function(numerator=formula, denominator, design, separate=FALSE,na.rm=FALSE,formula,covmat=FALSE,deff=FALSE,...){

    if (separate){
      strats<-sort(unique(design$strata[,1]))
      if (!design$has.strata)
        warning("Separate and combined ratio estimators are the same for unstratified designs")
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$strata[,1] %in% s),
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
        numerator<-model.frame(numerator,design$variables,na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,design$variables,na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    numerator<-as.matrix(numerator)
    denominator<-as.matrix(denominator)
    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    nas<-!complete.cases(all)
    if ((na.rm==TRUE) && any(nas)){
      design<-design[!nas,]
      if (NROW(design$cluster) == NROW(all)){
        ## subset by zero weights
        all[nas,]<-1
        numerator[nas,]<-0
        denominator[nas,]<-1
      } else {
        ## subset by actually dropping rows
        all<-all[!nas,,drop=FALSE]
        numerator<-numerator[!nas,,drop=FALSE]
        denominator<-denominator[!nas,,drop=FALSE]
      }
    }
    allstats<-svytotal(all,design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))

    
    vars<-matrix(ncol=nd,nrow=nn)

    if (deff) deffs<-matrix(ncol=nd,nrow=nn)
    
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-svyrecvar(r*1/design$prob, design$cluster, design$strata, design$fpc,
                            postStrata=design$postStrata)
        if (deff){
          deffs[i,j]<-deff(svytotal(r,design,deff=TRUE))
        }
      }
    }
    if (covmat){
        ii<-rep(1:nn,nd)
        jj<-rep(1:nd,each=nn)
        allr<-sweep(numerator[,ii]-t(as.vector(rval$ratio)*t(denominator[,jj,drop=FALSE])),
                    2, colSums(denominator[,jj,drop=FALSE]/design$prob),"/")
        vcovmat<-svyrecvar(allr*1/design$prob, design$cluster, design$strata, design$fpc,
                           postStrata=design$postStrata)
        colnames(vcovmat)<-colnames(denominator)[ii]
        rval$vcov<-vcovmat
    }
    colnames(vars)<-colnames(denominator)
    rownames(vars)<-colnames(numerator)
    rval$var<-vars
    if (deff) attr(rval,"deff")<-deffs
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }
