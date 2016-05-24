uniah=function(formula, trt=NULL, data=NULL, shape='unimodal', mode='unknown', M=NULL, maxdec=3, maxiter=10^3, eps=10^-3){
  #1. load/check R packages
  
  #2. check input parameters
  #2.1. Null
  if( is.null(formula) )  stop("formular is requred")
  if( is.null(shape) )    stop("shape is required")
  if( is.null(mode) )     stop("mode is required")
  if( is.null(eps) )      stop("eps is required")
  if( is.null(maxiter) )  stop("maxiter is required")
  
  #2.2 Length
  if( length(shape)!=1 )    stop("shape must be scalar")
  if( length(mode)!=1 )     stop("mode must be scalar")
  if( length(eps)!=1 )      stop("eps must be scalar")
  if( length(maxiter)!=1 )  stop("maxiter must be scalar")
  
  #2.3 NA or Inf
  if( (is.na(eps)||is.infinite(eps)) )          stop("eps must be a finite value")
  if( (is.na(maxiter)||is.infinite(maxiter)) )  stop("maxiter must be a finite value")

  #2.4 class
  if(class(formula)!='formula') stop("formula argument is wrong")

  if(class(shape)!='character') stop("shape must be character")
  if(any(grep("uni",tolower(shape)))==TRUE){
    shape="unimodal"
  }else if(any(grep("ush",tolower(shape)))==TRUE){
    shape="ushape"
  }else if(any(grep("u-sh",tolower(shape)))==TRUE){
    shape="ushape"
  }
  if(shape!='unimodal' && shape!='ushape')
    stop("shape must be either unimodal or ushape")

  if(class(mode)!='character') stop("mode must be character")
  
  if(any(grep("unk",tolower(mode)))==TRUE){
    mode="unknown"
  }else if(any(grep("uk",tolower(mode)))==TRUE){
    mode="unknown"    
  }
  if(mode!='known' && mode!='unknown')
    stop("mode must be either known or unknown")

  if(mode=="known"){
    if( is.null(M) )                  stop("M is required for a known mode")  
    if( length(M)!=1 )                stop("M is required for a known mode. M must be scalar")  
    if( (is.na(M)||is.infinite(M)) )  stop("M is required for a known mode. M must be a finite value")  
    if( !(is.numeric(M)) )            stop("M is required for a known mode. M must be numeric")    
  }

  if( !(is.numeric(eps)) )     stop("eps must be numeric")
  if( !(is.numeric(maxiter)) ) stop("maxiter must be be numeric")
  
  #2.5 eps & maxiter
  maxiter=ceiling(maxiter)
  if(eps<0)       stop("eps must be greater than or equal to zero")
  if(maxiter<=0)  stop("maxiter must be greater than zero")
  
  #3. load data with formula 
  #3.1 Surv response and cov
  mf=model.frame(formula=formula, data=data) #data must be data.frame
  
  #3.2 outcome
  surv.y=model.response(mf)
  if( !(is.Surv(surv.y)) ) stop("Response must be a survival outcome")
  
  if(ncol(surv.y)==2){ #for time-indepnedent
    type='ti'
    TIME=surv.y[,1];     STATUS=surv.y[,2]
  }else if(ncol(surv.y)==3){ #for time-depnedent (or time-independent)
    type='td'
    START=surv.y[,1]; STOP=surv.y[,2];  STATUS=surv.y[,3]
  }
  
  #3.3 cov
  Z=model.matrix(formula, data = mf, rhs = 2)
  if (ncol(Z)>=3) stop("The current version of isoph does not support multivariate explanatory variables")
  Z=Z[,2] #Z[,1] is intercept;  
  
  #3.4. check data entry
  if(type=='ti'){
    #NA or Inf
    if ( any(is.na(TIME)+is.na(STATUS)+is.na(Z)) ) stop("Data included NA")
    if ( any(is.infinite(TIME)+is.infinite(STATUS)+is.infinite(Z)) )  stop("Data included infinite values")
    
    #length
    if( !(length(TIME)==length(STATUS) & length(STATUS)==length(Z)) ) stop("Lengths of data in the formula argument are not matched")
    
    #right censored data
    if( min(TIME)<=0 ) stop("Time must be greater than zero")
    
  }else if(type=='td'){
    #NA or Inf
    #if ( any(is.na(START)+is.na(STOP)+is.na(STATUS)+is.na(Z)) ) stop("Data included NA")
    #if ( any(is.infinite(START)+is.infinite(STOP)+is.infinite(STATUS)+is.infinite(Z)) ) stop("Data included infinite values") 
    
    #length
    #if( !(length(START)==length(STOP) & length(STOP)==length(STATUS) & length(STATUS)==length(Z)) ) stop("Lengths of data are not matched")
    
    #right censored data
    #if( min(STOP)<=0 ) stop("Time must be greater than zero")
  }
  
  #3.5 X
  if(!is.null(data))    trt=data$trt
  if(type=='ti'){
    if(!is.null(trt)){
      if( any(is.na(trt)) ) stop("trt included NA")
      if( any(is.infinite(trt)) ) stop("trt included NA")
      
      uniq.trt=unique(trt)
      if(length(uniq.trt)!=2)  stop("trt must be coded by 0 and 1")
      if(min(uniq.trt)!=0)       stop("trt must be coded by 0 and 1")
      if(max(uniq.trt)!=1)       stop("trt must be coded by 0 and 1")
    }
  }else if(type=='td'){
  }
  
  
  #Censoring
  if ( length(unique(STATUS))>=3 ) stop("status has to be either 0 or 1")
  if ( !all(STATUS %in% c(0,1)) )  stop("status has to be either 0 or 1")
  
  if(sum(STATUS)<=2) stop("At least more than two numbers of event are needed.")
  
  #4. uniah
  if(type=='ti'){
    if(mode=='known'){
      est=uniah.ti.known(TIME=TIME, STATUS=STATUS, Z=Z, X=trt, shape=shape, K=M, maxdec=maxdec, maxiter=maxiter, eps=eps)
    }else if(mode=='unknown'){
      est=uniah.ti.unknown(TIME=TIME, STATUS=STATUS, Z=Z, X=trt, shape=shape, maxdec=maxdec, maxiter=maxiter, eps=eps)
    }
  }else if(type=='td'){
    stop("interval data (ot time-dependent covariate) is not supported for the current version of the uniah function")
  }

  est$call=match.call()
  est$formula=formula
  
  class(est)="uniah"
  
  est
}
