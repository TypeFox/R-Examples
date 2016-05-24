
regcalibrate<-function(design, ...) UseMethod("regcalibrate")

is.calibrated<-function(design){ !is.null(design$postStrata)}

##
## unbounded linear calibration using qr decomposition: less sensitive to
## collinearity than Deville & Sarndal's Newton algorithm.
##
regcalibrate.survey.design2<-function(design, formula, population,
                                   stage=NULL,  lambda=NULL, aggregate.stage=NULL,...){
  
  if (is.null(stage))
    stage<-if (is.list(population)) 1 else 0

  if (!is.null(aggregate.stage)){
    aggindex<-design$cluster[[aggregate.stage]]
  }
  
  if(stage==0){
    ## calibration to population totals
    mm<-model.matrix(formula, model.frame(formula, model.frame(design)))
    ww<-weights(design)
    if (is.null(lambda))
      sigma2<-rep(1,nrow(mm))
    else
      sigma2<-drop(mm%*%lambda)

    if (!is.null(aggregate.stage)){
      mm<-apply(mm,2,function(mx) ave(mx,aggindex))
      ww<-ave(ww,aggindex)
      sigma2<-ave(sigma2,aggindex)
    }
    whalf<-sqrt(ww)
    sample.total<-colSums(mm*ww)

    if(any(sample.total==0)){
      ## drop columsn where all sample and population are zero
      zz<-(population==0) & (apply(mm,2,function(x) all(x==0)))
      mm<-mm[,!zz]
      population<-population[!zz]
      sample.total<-sample.total[!zz]
    }

    
    if (length(sample.total)!=length(population))
      stop("Population and sample totals are not the same length.")
    
    if (!is.null(names(population)) && any(names(sample.total)!=names(population)))
      warning("Sample and population totals have different names.")

    tqr<-qr(mm*whalf/sqrt(sigma2))

    ## not needed
    ##if (is.null(lambda) && !all(abs(qr.resid(tqr,whalf*sigma2)/sigma2) <1e-5))
    ##  warning("Calibration models with constant variance must have an intercept")

    g<-rep(1,NROW(mm))

    Tmat<-crossprod(mm*whalf/sqrt(sigma2))
    
    tT<-solve(Tmat,population-sample.total)
    
    g<-drop(1+mm%*%tT/sigma2)
    
    
    design$prob<-design$prob/g
    
    caldata<- list(qr=tqr, w=g*whalf*sqrt(sigma2), stage=0, index=NULL)
    
  } else {
    ## Calibration within clusters (Sarndal's Case C)
    if (stage>NCOL(design$cluster))
      stop("This design does not have stage",stage)

    if (!is.null(aggregate.stage)){
      stop("aggregate= not implemented for calibration within clusters")
    }

    if (!all(length(population[[1]])==sapply(population,length)))
      stop("Population totals are not all the same length")
    
    clusters<-unique(design$cluster[,stage])
    nc<-length(clusters)
    
    caldata<-list(qr=vector("list",nc), w=vector("list",nc),
                  stage=stage,index=as.character(clusters))

    mm<-model.matrix(formula, model.frame(formula, model.frame(design)))

    if (is.null(lambda))
      sigma2<-rep(1,nrow(mm))
    else
      sigma2<-drop(mm%*%lambda)
    
    if(NCOL(mm)!=length(population[[1]]))
        stop("Population and sample totals are not the same length.")
      
    if (any(colnames(mm)!=names(population[[1]])))
      warning("Sample and population totals have different names.")
  
    stageweights<-1/apply(design$allprob[,1:stage,drop=FALSE],1,prod)
    if (any(duplicated(design$cluster[!duplicated(stageweights),stage])))
      stop("Weights at stage", stage, "vary within sampling units")

    cwhalf<-sqrt(weights(design)/stageweights)
    dwhalf<-sqrt(weights(design))
    tqr<-qr(mm)

    ## not needed
    ## if (is.null(lambda) && !all(abs(qr.resid(tqr,sigma2)) <1e-3))
    ##  stop("Calibration models with constant variance must have an intercept")
 
    for (i in 1:length(clusters)){ 
      cluster<-clusters[[i]]
      these<-which(cluster ==  as.character(design$cluster[,stage]))
      mmi<-mm[these,,drop=FALSE]
      sample.total<-colSums(mmi*cwhalf[these]*cwhalf[these])
      
      if(any(sample.total==0)){
        ## drop columsn where all sample and population are zero
        zz<-(population[[i]]==0) & (apply(mmi,2,function(x) all(x==0)))
        mmi<-mmi[,!zz,drop=FALSE]
        population[[i]]<-population[[i]][!zz]
        sample.total<-sample.total[!zz]
      }

      tqr<-qr(mmi*cwhalf[these]/sqrt(sigma2[these]))
      Tmat<-crossprod(mmi*cwhalf[these]/sqrt(sigma2[these]))
      tT<-solve(Tmat,population[[i]]-sample.total)
      g<-drop(1+mmi%*%tT/sigma2[these])
      design$prob[these]<-design$prob[these]/g
      caldata$qr[[i]]<-tqr
      caldata$w[[i]]<-g*stageweights[these]*sqrt(sigma2[these])*cwhalf[these]^2
    }
  }  
  class(caldata)<-"greg_calibration"
  
  design$postStrata<-c(design$postStrata, list(caldata))
  design$call<-sys.call(-1)
  
  design
}


regcalibrate.svyrep.design<-function(design, formula, population,compress=NA,lambda=NULL,
                                  aggregate.index=NULL,...){
  mf<-model.frame(formula, design$variables)
  mm<-model.matrix(formula, mf)
  ww<-design$pweights
  if (is.null(lambda))
    sigma2<-rep(1,nrow(mm))
  else
    sigma2<-drop(mm%*%lambda)
  
  repwt<-as.matrix(design$repweights)
  if (!design$combined.weights)
    repwt<-repwt*design$pweights

  if (inherits(aggregate.index,"formula")){
    if (length(aggregate.index)!=2)
      stop("aggregate.index must be a one-sided formula")
    aggregate.index<-model.frame(aggregate.index, design$variables)
    if (NCOL(aggregate.index)>1)
      stop("aggregate.index must specify a single variable")
    aggregate.index<-aggregate.index[[1]]
  }
  
  if (!is.null(aggregate.index)){
    if (sqrt(max(ave(ww,aggregate.index,FUN=var),na.rm=TRUE))>1e-2*mean(ww))
      warning("Sampling weights are not constant within clusters defined by aggregate.index")
    mm<-apply(mm,2,function(mx) ave(mx,aggregate.index))
    ww<-ave(ww,aggregate.index)
    sigma2<-ave(sigma2,aggregate.index)
    repwt<-apply(repwt,2,function(wx) ave(wx, aggregate.index))
  }
  whalf<-sqrt(ww)
  
  sample.total<-colSums(mm*ww)
  
  if(any(sample.total==0)){
    ## drop columsn where all sample and population are zero
    zz<-(population==0) & (apply(mm,2,function(x) all(x==0)))
    mm<-mm[,!zz]
    population<-population[!zz]
    sample.total<-sample.total[!zz]
  }
  
  if (length(sample.total)!=length(population))
    stop("Population and sample totals are not the same length.")
  if (!is.null(names(population)) && any(names(sample.total)!=names(population)))
    warning("Sample and population totals have different names.")
  
  Tmat<-crossprod(mm*whalf/sqrt(sigma2))
  
  tT<-solve(Tmat,population-sample.total)
  
  gtotal<-drop(1+mm%*%tT/sigma2)
  design$pweights<-design$pweights*gtotal
  
  for(i in 1:NCOL(repwt)){
    whalf<-sqrt(repwt[,i])
    Tmat<-crossprod(mm*whalf/sqrt(sigma2))
    sample.total<-colSums(mm*whalf*whalf)
    g<-drop(1+mm%*%solve(Tmat,population-sample.total)/sigma2)
    repwt[,i]<-as.vector(design$repweights[,i])*g
  }

  if (!design$combined.weights)
      repwt<-repwt/gtotal

  if (compress ||
      (is.na(compress && inherits(design$repweights,"repweights_compressed")))){
    repwt<-compressWeights(repwt)
  }
    
  design$repweights<-repwt
  design$call<-sys.call(-1)
  design$degf<-NULL
  design$degf<-degf(design)

  design
}
