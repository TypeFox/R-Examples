# In this function we consider different kind estimation procedures for COGARCH(P,Q)
# Model. 
is.COGARCH <- function(obj){
  if(is(obj,"yuima"))
    return(is(obj@model, "yuima.cogarch"))
  if(is(obj,"yuima.cogarch"))
    return(is(obj, "yuima.cogarch"))
  return(FALSE)
}
yuima.acf<-function(data,lag.max, forward=TRUE){
  if(forward==FALSE){  
    burndata<-lag.max+1
    dataused<-as.list(numeric(length=burndata+1))
    Time<-length(data)
    dataformean<-data[burndata:Time]
    dataused[[1]]<-mean(dataformean)
    dataused[[2]]<-mean(dataformean^2)
    mu<-mean(dataformean)
    leng <-length(dataformean)
    
    var1<-sum((dataformean-mu)*(dataformean-mu))/leng
    
    
    res1<-numeric(length=lag.max)
    elem<-matrix(0,lag.max,(Time-lag.max))
    for (t in (lag.max+1):(Time)){
  
  #     h<-leng-lag.max
  #     elem<-(data[(1+t):leng]-mu)*(data[1:h]-mu)/(leng*var1) 
  #     res1[t+1]<-sum(elem) 
      h <-c(lag.max:1)
      elem[,(t-lag.max)]<-(data[(t-h)[order(t-h,decreasing=TRUE)]]-mu)*(data[t]-mu)/(var1)
    }
    for(h in 3:(lag.max+2)){
      dataused[[h]]<-sum(elem[h-2,])/leng
    }  
    
    elem0<-rbind(t(as.matrix(dataformean)),elem)
  #   res1<-res1
  #   acfr<-res1[2:(lag.max+1)]  #analogously to Matlab program
  }else{
    burndata<-lag.max
    dataused<-as.list(numeric(length=burndata+1))
    Time<-length(data)
    dataformean<-data[1:(Time-burndata)]
    dataused[[1]]<-mean(dataformean)
    dataused[[2]]<-mean(dataformean^2)
    mu<-mean(dataformean)
    leng <-length(dataformean)
    
    var1<-sum((dataformean-mu)*(dataformean-mu))/leng
    
    
    res1<-numeric(length=lag.max)
    elem<-matrix(0,lag.max,(Time-lag.max))
    for (t in 1:(Time-(lag.max))){
      
      #     h<-leng-lag.max
      #     elem<-(data[(1+t):leng]-mu)*(data[1:h]-mu)/(leng*var1) 
      #     res1[t+1]<-sum(elem) 
      h <-c(1:lag.max)
      elem[,t]<-(data[(t+h)]-mu)*(data[t]-mu)/(var1)
    }
    for(h in 3:(lag.max+2)){
      dataused[[h]]<-sum(elem[h-2,])/leng
    }  
    
    elem0<-rbind(t(as.matrix(dataformean)),elem)
  }
  return(list(dataused=dataused, elem=elem0, leng=leng))
}


# The estimation procedure for cogarch(p,q) implemented in this code are based on the 
# Chadraa phd's thesis
gmm<-function(yuima, data = NULL, start, method="BFGS", fixed = list(), 
                           lower, upper, lag.max = NULL, equally.spaced = FALSE, aggregation=TRUE,
                           Est.Incr = "NoIncr", objFun = "L2"){
  print <- FALSE

  aggr.G <- equally.spaced
  call <- match.call()
  
  if(objFun=="L1" && method!="Nelder-Mead"){
    yuima.warn("Mean absolute error minimization is available only for 'method=Nelder-Mead'. yuima sets automatically 'method=Nelder-Mead'  ")
    method<-"Nelder.Mead"
  }
  
  if(objFun=="L1" && (length(fixed)!=0 || !missing(lower) || !missing(upper))){
    yuima.stop("Constraints are not allow for the minimization of Mean absolute error")
  }
  
  codelist.objFun <- c("L1","L2","L2CUE", "TWOSTEPS")
  if(any(is.na(match(objFun,codelist.objFun)))){
    yuima.stop("Value of objFun not available. Please choose among L1, L2, L2CUE, TWOSTEPS")
  }
   
  codelist.Est.Incr <- c("NoIncr","Incr","IncrPar")
  if(any(is.na(match(Est.Incr,codelist.Est.Incr)))){
    yuima.stop("Value of Est.Incr not available. Please choose among NoIncr, Incr, IncrPar ")
  }
  if( missing(yuima))
    yuima.stop("yuima object is missing.")
  
  if( missing(start) ) 
    yuima.stop("Starting values for the parameters are missing.")
  
  if( !is.COGARCH(yuima) )
    yuima.stop("The model is not an object of class yuima.")

  if( !is(yuima,"yuima") && missing(data) )
    yuima.stop("data are missing.")
  
  
  if(is(yuima,"yuima")){
    model<-yuima@model
    if(is.null(data)){
      observ<-yuima@data
    }else{observ<-data}
  }else{
    if(is(yuima,"yuima.cogarch")){
      model<-yuima
      if(is.null(data)){
        yuima.stop("Missing data")
      }
      observ<-data
    }
  }
  
  if(!is(observ,"yuima.data")){
    yuima.stop("Data must be an object of class yuima.data-class")  
  } 
  
  if( !missing(upper) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The upper requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }

  if( !missing(lower) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The lower constraints requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }

  if( !missing(fixed) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The fixed constraints requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }
  
  # We identify the model parameters
  info <- model@info
  
  numb.ar <- info@q
  ar.name <- paste(info@ar.par,c(numb.ar:1),sep="")
  
  numb.ma <- info@p
  ma.name <- paste(info@ma.par,c(1:numb.ma),sep="")
  
  loc.par <- info@loc.par
  #xinit.name <- paste(info@Latent.var, c(1:numb.ar), sep = "")

  xinit.name0 <- model@parameter@xinit
  idx <- na.omit(match(c(loc.par, ma.name), xinit.name0))
  xinit.name <- xinit.name0[-idx]

  meas.par <- model@parameter@measure
  
  if(length(meas.par)==0 && Est.Incr=="IncrPar"){
    yuima.warn("The dimension of measure parameters is zero, yuima changes 'Est.Incr = IncrPar' into 'Est.Incr = Incr'")
    Est.Incr <- "Incr"
  }
    
    
    
  fixed.name <- names(fixed)
  if(info@q==1){
  #  nm <- c(names(start), "EL1", "phi1","phi2")
    nm <- c(names(start), "EL1")
  }else{
  #  nm <- c(names(start), "EL1", "M2Lev","M4Lev")
    nm <- c(names(start), "EL1")
  }
  # We identify the index of parameters
  if(length(meas.par)!=0){
    if(info@q==1){  
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, measure.par, "EL1", "phi1","phi2")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, meas.par, "EL1")
    }else{
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, measure.par, "EL1", "M2Lev","M4Lev")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, meas.par, "EL1")
    }
  }else{
    if(info@q==1){
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1", "phi1","phi2")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1")
    }else{
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1", "M2Lev","M4Lev")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1")
    }
  }
  oo <- match(nm, fullcoeff)

  if(any(is.na(oo)))
    yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
  if(info@q==1){
  #  start <- c(start, EL1 = 1, phi1=-1, phi2=-1)
    start <- c(start, EL1 = 1)
  }else{
  #  start <- c(start, EL1 = 1, M2Lev=1, M4Lev=2)
    start <- c(start, EL1 = 1)
  }
  start <- start[order(oo)]
  
  nm <- names(start)

  ar.idx <- match(ar.name, fullcoeff)
  ma.idx <- match(ma.name, fullcoeff)
  loc.idx <- match(loc.par, fullcoeff)
  meas.idx <- match(meas.par, fullcoeff)
  fixed.idx <- match(fixed.name, fullcoeff)
  EL1.idx <- match("EL1",fullcoeff)  # We decide to pass EL1 as parameter !!!

  env <- new.env()
#  n <- length(observ)[1]

  #n <- attr(observ@original.data,"tsp")[2]
  
  n <- length(index(observ@original.data))

  #Lag
  assign("lag", lag.max, envir=env)
  
  # Data
  assign("Data",  as.matrix(onezoo(observ)[,1]), envir=env)
  assign("deltaData",  n/index(observ@zoo.data[[1]])[n], envir=env)
  assign("time.obs",length(env$Data),envir=env)
  
  
  # Order
  assign("p", info@p, envir=env)
  assign("q", info@q, envir=env)

  # Idx
  assign("ar.idx", ar.idx, envir=env)
  assign("ma.idx", ma.idx, envir=env)
  assign("loc.idx", loc.idx, envir=env)
  assign("meas.idx", meas.idx, envir=env)
  assign("EL1.idx", EL1.idx, envir=env)
  
  objFunDummy <- NULL
  if(objFun=="TWOSTEPS"||objFun=="L2CUE"){
    objFunDummy <- objFun
    objFun <- "L2"
    
  }
  assign("objFun",objFun, envir=env)
  
  if(aggr.G==TRUE){
    if(floor(n/index(observ@zoo.data[[1]])[n])!=env$deltaData){
      yuima.stop("the n/Terminal in sampling information is not an integer. equally.spaced=FALSE is recommended")
    }
  }
  
  if(aggr.G==TRUE){
    # Aggregate returns G
    #dt<-round(deltat(onezoo(observ)[,1])*10^5)/10^5
    # Time<-index(observ@zoo.data[[1]])[n]
    G_i <- diff(env$Data[seq(1,length(env$Data),by=env$deltaData)])
    r<-1
  }else{
    dummydata<-index(onezoo(observ)[,1])
    unitarytime<-floor(dummydata)
    index<-!duplicated(unitarytime)
    G_i <- diff(env$Data[index]) 
    
    r <- 1
  }
  d <- min(floor(sqrt(length(G_i))),env$lag)
  assign("d", d, envir=env)
  typeacf <- "correlation"
  assign("typeacf", typeacf, envir=env)
# CovQuad <- acf(G_i^2,plot=FALSE,lag.max=d,type=typeacf)$acf[-1]
  
  example<-yuima.acf(data=G_i^2,lag.max=d)
  dummyEmpiricalMoM<-as.numeric(example$dataused)
  CovQuad <- as.numeric(example$dataused)[-c(1:2)]
#   
  
  
  
  assign("G_i", G_i, envir=env)
  assign("r", r, envir=env)
  #mu_G2 <- as.numeric(example$dataused)[1]
  mu_G2<-mean(G_i^2)
  assign("mu_G2", mu_G2, envir=env)
  #var_G2 <- as.numeric(example$dataused)[2] - mu_G2^2
  var_G2 <- mean(G_i^4) - mu_G2^2
  assign("var_G2", var_G2, envir=env)
  assign("score",example$elem,envir=env )
  assign("leng",example$leng,envir=env )
  #CovQuad <-log(abs(yuima.acf(data=G_i^2,lag.max=min(d,env$lag))))
  assign("CovQuad", CovQuad, envir=env)

objectiveFun <- function(p,env) {
  mycoef <- as.list(p)
    if(length(c(fixed.idx, meas.idx))>0){ ## SMI 2/9/14
      names(mycoef) <- nm[-c(fixed.idx,meas.idx)] ## SMI 2/9/14
    }else{
      names(mycoef) <- nm
   } 
  ErrTerm(yuima=yuima, param=mycoef, print=print, env)
}
if(method!="L-BFGS-B"&&method!="brent"){
  out<- optim(start, objectiveFun, method = method, env=env, hessian = TRUE)
}else{
  if(length(fixed)!=0 && !missing(upper) && !missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                fixed=as.numeric(fixed), 
                lower=as.numeric(lower), 
                upper=as.numeric(upper), env=env)
  }else{
    if(!missing(upper) && !missing(lower)){
      out<- optim(start, objectiveFun, method = method, 
                  lower=as.numeric(lower), 
                  upper=as.numeric(upper), env=env)
    }
    if(length(fixed)!=0 && !missing(lower)){
      out<- optim(start, objectiveFun, method = method, 
                  fixed=as.numeric(fixed), 
                  lower=as.numeric(lower), env=env)
    }
    if(!missing(upper) && length(fixed)!=0){
      out<- optim(start, objectiveFun, method = method, 
                  fixed=as.numeric(fixed), 
                  upper=as.numeric(upper), env=env)
    }
  }
  
  if(length(fixed)!=0 && missing(upper) && missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                fixed=as.numeric(fixed), env=env)
  }
  
  if(length(fixed)==0 && !missing(upper) && missing(lower)){
    out<- optim(start, objectiveFun, method = method,
                upper=as.numeric(upper), env=env)
  }
  
  if(length(fixed)==0 && missing(upper) && !missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                lower=as.numeric(lower), env=env)
  }
  
}

# Alternative way for calculating Variance Covariance Matrix

# sig2eps<-out$value/d
# 
# dumHess<- out$hessian[c(ar.name, ma.name),c(ar.name, ma.name)]/(2*sig2eps)
# vcovgmm<-solve(dumHess)
# sqrt(diag(vcovgmm))
 if(!is.null(objFunDummy)){
   assign("objFun",objFunDummy, envir=env)
   
   bvect<-out$par[ar.name]
   bq<-bvect[1]
   avect<-out$par[ma.name]
   a1<-avect[1]
   
   out$par[loc.par]<-(bq-a1)*mu_G2/(bq*r)
   
   # Determine the Variance-Covariance Matrix
   if(length(meas.par)!=0){
     idx.dumm<-match(meas.par,names(out$par))
     out$par<-out$par[- idx.dumm]
   }
   dimOutp<-length(out$par)-(1+info@q)
   coef <- out$par[c(1:dimOutp)] 
   vcov<-matrix(NA, dimOutp, dimOutp)
   names_coef<-names(coef)
   colnames(vcov) <- names_coef
   rownames(vcov) <- names_coef
   mycoef <- start
   min <- out$value 
   # # call
     gradVect0<-MM_grad_Cogarch(p=info@p, q=info@q, 
                                acoeff=avect,cost=out$par[loc.par], b=bvect,  
                                r=env$r, h=seq(1, env$d, by = 1)*env$r, type=typeacf, 
                                m2=env$mu_G2, var=env$var_G2)
     
     score0 <- MM_Cogarch(p=info@p, q=info@q,
                          acoeff=avect,cost=out$par[loc.par], b=bvect,  
                          r=env$r, h=seq(1, env$d, by = 1)*env$r, type=typeacf, 
                          m2=env$mu_G2, var=env$var_G2)

     idx.aaa<-match(loc.par,names_coef)           
     #  gradVect <- gradVect0[names_coef[-idx.aaa], ]
     gradVect <- gradVect0[names_coef[-idx.aaa],CovQuad>0]
     #  score <- c(score0$acfG2)%*%matrix(1,1,example$leng)
     score <- c(score0$acfG2[CovQuad>0])%*%matrix(1,1,example$leng)
     
     #We need to write the matrix W for the matrix sandwhich
     #plot(as.numeric(example$dataused)[-1],type="h")
     #S_matrix
     
     #  EmpirScore <-score-example$elem[-1,]
     exampelem <-example$elem[-1,]
     EmpirScore <-score-exampelem[CovQuad>0,]
     
     Omega_est <- tryCatch((1/example$leng*EmpirScore%*%t(EmpirScore)), 
                         error=function(theta){NULL})
     W_est <-tryCatch(chol2inv(Omega_est),error=function(theta){NULL})
   if(!is.null(W_est)){
     assign("W_est",W_est,envir=env)
     start0<-unlist(start)
     start0[names(out$par)]<-out$par
     start <- as.list(start0)
     #start<-out$par
    if(method!="L-BFGS-B"&&method!="brent"){
     out<- optim(start, objectiveFun, method = method, env=env)
   }else{
     if(length(fixed)!=0 && !missing(upper) && !missing(lower)){
       out<- optim(start, objectiveFun, method = method, 
                   fixed=as.numeric(fixed), 
                   lower=as.numeric(lower), 
                   upper=as.numeric(upper), env=env)
     }else{
       if(!missing(upper) && !missing(lower)){
         out<- optim(start, objectiveFun, method = method, 
                     lower=as.numeric(lower), 
                     upper=as.numeric(upper), env=env)
       }
       if(length(fixed)!=0 && !missing(lower)){
         out<- optim(start, objectiveFun, method = method, 
                     fixed=as.numeric(fixed), 
                     lower=as.numeric(lower), env=env)
       }
       if(!missing(upper) && length(fixed)!=0){
         out<- optim(start, objectiveFun, method = method, 
                     fixed=as.numeric(fixed), 
                     upper=as.numeric(upper), env=env)
       }
     }
     
     if(length(fixed)!=0 && missing(upper) && missing(lower)){
       out<- optim(start, objectiveFun, method = method, 
                   fixed=as.numeric(fixed), env=env)
     }
     
     if(length(fixed)==0 && !missing(upper) && missing(lower)){
       out<- optim(start, objectiveFun, method = method,
                   upper=as.numeric(upper), env=env)
     }
     
     if(length(fixed)==0 && missing(upper) && !missing(lower)){
       out<- optim(start, objectiveFun, method = method, 
                   lower=as.numeric(lower), env=env)
     }
     
   }
   }else{
     yuima.warn("Method TWOSTEPS or L2CUE Changed in L2 since First W failed to compute")
   }
 }
 bvect<-out$par[ar.name]
 bq<-bvect[1]
 avect<-out$par[ma.name]
 a1<-avect[1]
 
  out$par[loc.par]<-(bq-a1)*mu_G2/(bq*r)

 # Determine the Variance-Covariance Matrix
 if(length(meas.par)!=0){
   idx.dumm<-match(meas.par,names(out$par))
   out$par<-out$par[- idx.dumm]
 }
 dimOutp<-length(out$par)-(1+info@q)
 coef <- out$par[c(1:dimOutp)] 
 vcov<-matrix(NA, dimOutp, dimOutp)
  names_coef<-names(coef)
  colnames(vcov) <- names_coef
  rownames(vcov) <- names_coef
 # mycoef <- start
  min <- out$value 
# # call
  if(objFun!="L1"){
        gradVect0<-MM_grad_Cogarch(p=info@p, q=info@q, 
                  acoeff=avect,cost=out$par[loc.par], b=bvect,  
                  r=env$r, h=seq(1, env$d, by = 1)*env$r, type=typeacf, 
                  m2=env$mu_G2, var=env$var_G2)
  
        score0 <- MM_Cogarch(p=info@p, q=info@q,
                         acoeff=avect,cost=out$par[loc.par], b=bvect,  
                         r=env$r, h=seq(1, env$d, by = 1)*env$r, type=typeacf, 
                         m2=env$mu_G2, var=env$var_G2)
    if(objFun == "L2"){
    #  min <- log(sum((score0$acfG2[CovQuad>0]-CovQuad[CovQuad>0])^2))
      
      min <- log(sum((score0$acfG2-CovQuad)^2))
      #min <- log(sum((score0$acfG2[CovQuad>0]-CovQuad[CovQuad>0])^2))
    }          
  idx.aaa<-match(loc.par,names_coef)           
  #  gradVect <- gradVect0[names_coef[-idx.aaa], ]
   gradVect <- gradVect0[names_coef[-idx.aaa],CovQuad>0]
  #  score <- c(score0$acfG2)%*%matrix(1,1,example$leng)
  score <- c(score0$acfG2[CovQuad>0])%*%matrix(1,1,example$leng)
  
  #We need to write the matrix W for the matrix sandwhich
#plot(as.numeric(example$dataused)[-1],type="h")
  #S_matrix

#  EmpirScore <-score-example$elem[-1,]
    exampelem <-example$elem[-1,]
    EmpirScore <-score-exampelem[CovQuad>0,]
  
  Omega_est<-tryCatch((1/example$leng*EmpirScore%*%t(EmpirScore)), 
                   error=function(theta){NULL})
  if(is.null(Omega_est)){
    Omega_est<-matrix(NA,dim(EmpirScore)[1],dim(EmpirScore)[1])
  }
 if(is.null(objFunDummy)){
    Gmatr<-gradVect%*%t(gradVect)
    CentMat<-gradVect%*%Omega_est%*%t(gradVect)
    Var_Matr0 <- tryCatch(solve(Gmatr)%*%CentMat%*%solve(Gmatr)/example$leng,
                          error=function(theta){NULL})
 }else{
   #Gmatr<-gradVect%*%W_est%*%t(gradVect)
   #CentMat<-gradVect%*%W_est%*%Omega_est%*%W_est%*%t(gradVect)
   Var_Matr0 <- tryCatch(solve(gradVect%*%solve(Omega_est)%*%t(gradVect))/example$leng,
                         error=function(theta){NULL})
   #Var_Matr0 <- solve(Gmatr)/example$leng
 }

  
  if(!is.null(Var_Matr0)){
    aaa<-dimOutp-1
    vcov[c(1:aaa),c(1:aaa)]<-Var_Matr0  
  }

#    Var_Matr <- solve(gradVect%*%t(gradVect))
  #vcov <- Var_Matr
#    vcov[loc.par,]<-NA
#    vcov[,loc.par]<-NA
#out$par[loc.par]<-(bq-a1)*mu_G2/(bq*r)
  }

  
  # Build an object of class mle
  if(Est.Incr=="NoIncr"){
      res<-new("cogarch.gmm", call = call, coef = coef, fullcoef = unlist(coef), 
                vcov = vcov, min = min, details = list(), 
                method = character(),
                model = model,
                objFun = objFun 
               )
  }
  if(Est.Incr=="Incr"||Est.Incr=="IncrPar"){
    L.Incr<-cogarchNoise(yuima.cogarch=model, data=observ, 
                            param=as.list(coef), mu=1)
    ttt<-observ@zoo.data[[1]]
    tt<-index(ttt)
    L.Incr_Fin <- zoo(L.Incr,tt[(1+length(tt)-length(L.Incr)):length(tt)])
  }
  if(Est.Incr=="Incr"){
  # Build an object of class cogarch.gmm.incr
      res<-new("cogarch.gmm.incr", call = call, coef = coef, fullcoef = unlist(coef), 
                vcov = vcov, min = min, details = list(), 
                method = character(),
                Incr.Lev = L.Incr_Fin,
                model = model, nobs=as.integer(length(L.Incr)+1),
                logL.Incr = numeric(),
                objFun= objFun
                 )
  }
  if(Est.Incr=="IncrPar"){
    #estimationLevy

    fixedCon <- constdum(fixed, meas.par)
    lowerCon <- constdum(lower, meas.par)
    upperCon <- constdum(upper, meas.par)
    if(aggregation==TRUE){
      if(floor(n/index(observ@zoo.data[[1]])[n])!=env$deltaData){
        yuima.stop("the n/Terminal in sampling information is not an integer. Aggregation=FALSE is recommended")
      }
      inc.levy1<-diff(cumsum(c(0,L.Incr))[seq(from=1,
                                                to=yuima@sampling@n[1],
                                                by=env$deltaData
      )])
    }else{
      inc.levy1 <- L.Incr
    }
    
    result.Lev <- gmm.Est.Lev(Increment.lev=c(0,inc.levy1), 
                                  param0=start[meas.par],
                                  fixed = fixedCon[meas.par],
                                  lower=lowerCon[meas.par],
                                  upper=upperCon[meas.par],
                                  measure=model@measure,
                                  measure.type=model@measure.type,
                                  aggregation=aggregation,
                                  dt=1/env$deltaData
                              )   
    
    if(is.null(result.Lev)){
       res<-new("cogarch.gmm.incr", call = call, coef = coef, fullcoef = unlist(coef), 
             vcov = vcov, min = min, details = list(), 
             method = character(),
             Incr.Lev=L.Incr_Fin,
             model = model, nobs=as.integer(length(L.Incr)+1),
             logL.Incr = numeric(),
             objFun= objFun
             )
    
    }
    else{
      Inc.Parm<-result.Lev$estLevpar
      IncVCOV<-result.Lev$covLev
      if(length(meas.par)==length(Inc.Parm)){
        names(Inc.Parm)<-meas.par
        rownames(IncVCOV)<-as.character(meas.par)
        colnames(IncVCOV)<-as.character(meas.par)
      }
      name.parm.cog<-names(coef)
      coef<-c(coef,Inc.Parm)
      
      names.par<-names(coef)
      cov<-NULL
      cov<-matrix(NA,length(names.par),length(names.par))
      rownames(cov)<-names.par
      colnames(cov)<-names.par
      
      cov[unique(name.parm.cog),unique(name.parm.cog)]<-vcov
      cov[names(Inc.Parm),names(Inc.Parm)]<-IncVCOV
      cov<-cov
      
      res<-new("cogarch.gmm.incr", call = call, coef = coef, fullcoef = unlist(coef), 
               vcov = cov, min = min, details = list(), 
               method = character(),
               Incr.Lev=L.Incr_Fin,
               model = model, nobs=as.integer(length(L.Incr)+1),
               logL.Incr = tryCatch(-result.Lev$value,error=function(theta){NULL}),
               objFun= objFun
               )
      
    }


 }
 return(res)
  
}

constdum<-function(fixed, meas.par){
  fixedCon<-list()
  if(!missing(fixed)){
    measure.par.dum.idx<-na.omit(names(fixed[meas.par]))
    if(length(measure.par.dum.idx)!=0){
      fixedCon<-fixed[measure.par.dum.idx]
    }
  }
}


gmm.Est.Lev<-function(Increment.lev, 
            param0,
            fixed=list(),
            lower=list(),
            upper=list(),
            measure,
            measure.type,
            aggregation,
            dt){
  
  fixed.carma <- unlist(fixed)    
  lower.carma <- unlist(lower)
  upper.carma <- unlist(upper)
  
  Dummy <- TRUE 
  
  
  CPlist <- c("dnorm","dgamma", "dexp")
  codelist <- c("rngamma","rNIG","rIG", "rgamma")
  if(measure.type=="CP"){
    tmp <- regexpr("\\(", measure$df$exp)[1]
    measurefunc <- substring(measure$df$exp, 1, tmp-1)
    if(is.na(match(measurefunc,CPlist))){
      yuima.warn("COGARCH(p,q): Other compound poisson processes will be implemented as soon as")
      Dummy <- FALSE 
    }
    
  }
  
  
  if(measure.type=="code"){
    tmp <- regexpr("\\(", measure$df$exp)[1]
    measurefunc <- substring(measure$df$exp, 1, tmp-1)
    if(is.na(match(measurefunc,codelist))){
      yuima.warn("COGARCH(p,q): Other Levy measure will be implemented as soon as possible")
      Dummy <- FALSE
      
    }      
  }  
  if(Dummy){
    #env$h<-dt
    result.Lev <- yuima.Estimation.Lev(Increment.lev=Increment.lev,
                                     param0=unlist(param0),
                                     fixed.carma=fixed.carma,
                                     lower.carma=lower.carma,
                                     upper.carma=upper.carma,
                                     measure=measurefunc,
                                     measure.type=measure.type,
                                     aggregation=aggregation,
                                     dt=dt)
  }else{
    result.Lev<-NULL
  }
  return(result.Lev)
}



ErrTerm <- function(yuima, param, print, env){
  typeacf <- env$typeacf
  param <- as.numeric(param)
  
  G_i <- env$G_i
  r <- env$r
  mu_G2 <- env$mu_G2
  var_G2 <- env$var_G2
  d <- env$d
  CovQuad <-env$CovQuad

  h <- seq(1, d, by = 1)*r
  cost <- env$loc.idx
  b <- env$ar.idx
  a <- env$ma.idx
  meanL1 <- param[env$EL1.idx]
#   meanL1 <- 1
#   if(env$q == 1){
# 
#   beta <- param[cost]*param[b]
#   eta <- param[b]
#   phi <- param[a]
#    
# #   phi1 <- param[env$phi1.idx]
# #   phi2 <- param[env$phi2.idx]
#   #theo_mu_G2 <- meanL1*r*beta/abs(phi1)
#   phi1 <- meanL1*r*beta/mu_G2
#   
#   termA <- (6*mu_G2/r*beta/abs(phi1)*(2*eta/phi-meanL1)*(r-(1-exp(-r*abs(phi1)))/abs(phi1))+2*beta^2/phi^2*r)
#   phi2 <-2*termA*abs(phi1)/((var_G2-2*mu_G2^2)*abs(phi1)+termA)
#   if(typeacf == "covariance"){
#   TheoCovQuad <- meanL1*beta^2/abs(phi1)^3*(2*eta/phi-meanL1)*
#      (2/abs(phi2)-1/abs(phi1))*(1-exp(-r*abs(phi1)))*(exp(r*abs(phi1))-1)*
#      exp(-h*abs(phi1))
#   }else{
#     TheoCovQuad <- meanL1*beta^2/abs(phi1)^3*(2*eta/phi-meanL1)*
#       (2/abs(phi2)-1/abs(phi1))*(1-exp(-r*abs(phi1)))*(exp(r*abs(phi1))-1)*
#       exp(-h*abs(phi1))/(var_G2)
#   }
#   
# } 
 if(env$q >= 1){
   TheoCovQuad <- numeric(length = length(h))
   cost<-param[cost]
#   for(i in c(1:length(h))){
      MomentCog <- MM_Cogarch(p = env$p, q = env$q,  acoeff=param[a],
                              cost=cost,
                              b=param[b],  r = r, h = h, 
                              type = typeacf, m2=mu_G2, var=var_G2)
   
   
   TheoCovQuad <- MomentCog$acfG2
 #  }
   theo_mu_G2 <- MomentCog$meanG2
   #param[cost]<-MomentCog$cost
 }

 if(env$objFun=="L2"){
#  res <- log(sum((TheoCovQuad[CovQuad>0]-CovQuad[CovQuad>0])^2))

#    emp <- log(CovQuad[CovQuad>0])
#    theo <- log(TheoCovQuad[CovQuad>0])
#   res <- log(sum((abs(TheoCovQuad)-abs(CovQuad))^2))
   
   
 res <- sum((log(TheoCovQuad[CovQuad>0])-log(CovQuad[CovQuad>0]))^2)

 #    res <- sum((TheoCovQuad[CovQuad>0]-CovQuad[CovQuad>0])^2)
 #  res <- sum((TheoCovQuad-CovQuad)^2)
  
#  res <- sum((log(abs(TheoCovQuad))-log(abs(CovQuad)))^2)
  
  # res <- sum((log(TheoCovQuad[CovQuad>0]))-log(CovQuad[CovQuad>0]))^2)
  return(res)
 }

if(env$objFun=="TWOSTEPS"){
  TheoMoM<-log(TheoCovQuad[CovQuad>0])
  #TheoMoM<-log(abs(TheoCovQuad))
  #TheoMoM<-TheoCovQuad[CovQuad>0]
  #dumScore <- (TheoMoM%*%matrix(1,1,env$leng))-env$score[-1,]
  #W_est<-chol2inv(dumScore%*%t(dumScore)/env$leng)
  W_est<-env$W_est
  CompMoM0<-TheoMoM-c(log(CovQuad[CovQuad>0]))
  #CompMoM0<-TheoMoM-c(log(abs(CovQuad)))
  #CompMoM0<-TheoMoM-CovQuad[CovQuad>0]
  CompMoM<-matrix(CompMoM0,1,length(CompMoM0))
  res <- as.numeric(CompMoM%*%W_est%*%t(CompMoM))
  
  
#   TheoMoM<-TheoCovQuad[CovQuad>0]
#   intrDum<-env$score[-1,]
#   
#   dumScore <- (TheoMoM%*%matrix(1,1,env$leng))-intrDum[CovQuad>0,]
#   W_est<-chol2inv(dumScore%*%t(dumScore)/env$leng)
#   CompMoM0<-TheoMoM-c(CovQuad[CovQuad>0])
#   CompMoM<-matrix(CompMoM0,1,length(CompMoM0))
#   res <- as.numeric(CompMoM%*%W_est%*%t(CompMoM))
  
  
  return(res)
}


  if(env$objFun=="L1"){
    res <- log(sum(abs(c(TheoCovQuad)-c(CovQuad))))
    return(res)
  }


 if(env$objFun=="L2CUE"){  
   #TheoMoM<-TheoCovQuad
   #
   #TheoMoM<-log(TheoCovQuad[CovQuad>0])
   TheoMoM<-TheoCovQuad[CovQuad>0]
   intrDum<-env$score[-1,]
   
   dumScore <- (TheoMoM%*%matrix(1,1,env$leng))-intrDum[CovQuad>0,]
   W_est<-chol2inv(dumScore%*%t(dumScore)/env$leng)
   CompMoM0<-TheoMoM-c(CovQuad[CovQuad>0])
   CompMoM<-matrix(CompMoM0,1,length(CompMoM0))
   res <- as.numeric(CompMoM%*%W_est%*%t(CompMoM))
   return(res)
 }

 
}

MM_Cogarch <- function(p, q, acoeff,cost, b,  r, h, type, m2, var){
  # The code developed here is based on the acf for squared returns derived in the 
  # Chaadra phd Thesis
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  a[1:p,1] <- acoeff

  bq <- b[1]
  a1 <- a[1]

# # Matching only the autocorrelation we are not able to estimate the a_0 parameter
# nevertheless using the theoretical formula of the expectaion of squared returns we
# are able to fix this parameter for having a perfect match  between theoretical and 
# empirical mean

# We recall that under the assumption of levy  process is centered  the mean at time 1 of the 
# squared returns and the mean of the corresponding levy measure are equals.

  mu<-1 # we assume variance of the underlying L\'evy is one
  meanL1<-mu
  

  cost<-(bq-mu*a1)*m2/(bq*r)
  B<- MatrixA(b[c(q:1)])
  B_tilde <- B+mu*e%*%t(a)
  meanG2 <- cost*bq*r/(bq-mu*a1)*meanL1
  Inf_eps <- IdentM <- diag(q)
  if(q==1){
    Inf_eps[q,q] <- -1/(2*B_tilde[1,1]) 
  }
  if(q==2){
    Inf_eps[q,q] <- -B_tilde[2,1]
    Inf_eps <- 1/(2*B_tilde[2,1]*B_tilde[2,2])*Inf_eps
  }
  if(q==3){
    Inf_eps[1,1] <- B_tilde[3,3]/B_tilde[3,1]
    Inf_eps[q,q] <- -B_tilde[3,2]
    Inf_eps[1,3] <- -1
    Inf_eps[3,1] <- -1
    Inf_eps <- 1/(2*(B_tilde[3,3]*B_tilde[3,2]+B_tilde[3,1]))*Inf_eps
  }
  if(q>=4){
    Inf_eps <- round(AsympVar(B_tilde,e,lower=0,upper=100,delta=0.01)*10^5)/10^5
  }
   
  #term <- expm(B_tilde*h)
  term  <- lapply(h,"yuimaExpm",B=B_tilde)
  #invB <- solve(B_tilde) # In this case we can use the analytical form for Companion Matrix???
  # We invert the B_tilde using the Inverse of companion matrix
  if(q>1){
  invB<-rbind(c(-B_tilde[q,-1],1)/B_tilde[q,1],cbind(diag(q-1),matrix(0,q-1,1)))
  }else{invB<-1/B_tilde}
  term1 <- invB%*%(IdentM-expm(-B_tilde*r))
  term2 <- (expm(B_tilde*r)-IdentM)
  
  P0_overRho <- 2*mu^2*(3*invB%*%(invB%*%term2-r*IdentM)-IdentM)%*%Inf_eps
  Q0_overRho <- 6*mu*((r*IdentM-invB%*%term2)%*%Inf_eps
                      -invB%*%(invB%*%term2-r*IdentM)%*%Inf_eps%*%t(B_tilde))%*%e
#   Q0_overRho <- 6*mu*((r*IdentM-invB%*%term2)%*%Inf_eps
#                     -invB%*%(invB%*%term2-r*IdentM)%*%Inf_eps%*%t(B))%*%e
  m_overRho <- as.numeric(t(a)%*%Inf_eps%*%a)
  Den<- (m_overRho*meanL1^2/m2^2*var*r^2+t(a)%*%Q0_overRho+t(a)%*%P0_overRho%*%a+1) 
  num <-(meanL1^2/m2^2*var-2*mu^2)*r^2
  rh0 <- as.numeric(num/Den)


  Inf_eps1 <- Inf_eps*rh0
  Ph_withouterm <- term1%*%invB%*%term2%*%Inf_eps1*mu^2
  Ph<-lapply(term,"%*%",Ph_withouterm)
  Qh_without <- term1%*%(-term2%*%Inf_eps1-invB%*%term2%*%Inf_eps1%*%t(B_tilde))%*%e*mu
  Qh<-lapply(term,"%*%",Qh_without) 
# Qh <- mu*term%*%term1%*%(-term2%*%Inf_eps1-invB%*%term2%*%Inf_eps1%*%t(B))%*%e
  m <- m_overRho*rh0
  atrasp<-t(a)
  if(type=="correlation"){
    VarTheo<-as.numeric(rh0*atrasp%*%P0_overRho%*%a+rh0*atrasp%*%Q0_overRho+2*r^2*mu^2+rh0)
    acfG2 <- (as.numeric(lapply(Ph,"yuimaQuadrForm",at=atrasp,a=a)) 
              + as.numeric(lapply(Qh,"yuimaInnerProd",at=atrasp)))/VarTheo
  }else{
    coeff <- cost^2*bq^2/((1-m)*(bq-mu*a1)^2)
    acfG2 <- coeff*( as.numeric(lapply(Ph,"yuimaQuadrForm",at=atrasp,a=a)) 
                     + as.numeric(lapply(Qh,"yuimaInnerProd",at=atrasp)) )
  }
  res <- list(acfG2=acfG2, meanG2=meanG2)
  return(res)
}


MM_grad_Cogarch <- function(p, q, acoeff,cost, b,  r, h, type, m2, var){
  eps<-10^(-3)
  PartialP<-matrix(0,p,length(h))
  epsA<-eps*diag(p)
  for(i in c(1:p)){
    LeftApproxP<-MM_Cogarch(p, q, acoeff-epsA[i,],cost, b,  r, h, type, m2, var)
    RightApproxP<-MM_Cogarch(p, q, acoeff[i]+epsA[i,],cost, b,  r, h, type, m2, var)
    PartialP[i,]<-(RightApproxP$acfG2-LeftApproxP$acfG2)/(2*eps)
  }
  PartialQ<-matrix(0,q,length(h))
  epsB<-eps*diag(q)
  for(i in c(1:q)){
    LeftApproxQ<-MM_Cogarch(p, q, acoeff,cost, b-epsB[i,],  r, h, type, m2, var)
    RightApproxQ<-MM_Cogarch(p, q, acoeff,cost, b+epsB[i,],  r, h, type, m2, var)
    PartialQ[i,]<-(RightApproxQ$acfG2-LeftApproxQ$acfG2)/(2*eps)
  }
#   LeftApproxCost<-MM_Cogarch(p, q, acoeff,cost-eps, b,  r, h, type, m2, var)
#   RightApproxCost<-MM_Cogarch(p, q, acoeff,cost+eps, b,  r, h, type, m2, var)
#   PartialCost0<-(c(RightApproxCost$meanG2,RightApproxCost$acfG2)-c(LeftApproxCost$meanG2,LeftApproxCost$acfG2))/(2*eps)
#   PartialCost<-matrix(PartialCost0,1, (length(h)+1))
  namesCoeff<-names(c(acoeff,b))
  res<-rbind(PartialP,PartialQ)
  rownames(res)<-namesCoeff
  return(res)
}

yuimaQuadrForm<-function(P,a,at){at%*%P%*%a}
yuimaInnerProd<-function(P,at){at%*%P}

yuimaExpm<-function(h,B){expm(B*h)}

AsympVar<-function(B_tilde,e,lower=0,upper=100,delta=0.1){
  part<-seq(lower,upper-delta, by=delta)+delta/2
  last <- length(part)
  Integrand <- as.matrix(matrix(0,length(e),length(e)))
  for(i in c(1:last)){
  Integrand<-Integrand+expm(B_tilde*part[i])%*%e%*%t(e)%*%expm(t(B_tilde)*part[i])*delta
  }
  return(Integrand)
}


setMethod("plot",signature(x="cogarch.gmm.incr"),
          function(x, type="l" ,...){
            Time<-index(x@Incr.Lev)
            Incr.L<-coredata(x@Incr.Lev)
            if(is.complex(Incr.L)){
              yuima.warn("Complex increments. We plot only the real part")
              Incr.L<-Re(Incr.L)
            }
#             model <- x@model
#             EndT <- Time[length(Time)]
#             numb <- (length(Incr.L)+1)
            plot(x= Time,y= Incr.L, type=type,...)
          }
)


setMethod("summary", "cogarch.gmm",
          function (object, ...)
          {
            cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
            obj<- object@min
            labFun <- object@objFun
            
            tmp <- new("summary.cogarch.gmm", call = object@call, coef = cmat,
                       m2logL = 0,
                       #model = object@model,
                       objFun = labFun,
                       objFunVal = obj
            )
            tmp
          }
)

# errorfun <- function(estimates, labelFun = "L2"){
#   if(LabelFun == "L2"){
#     
#   }
#     
# }

setMethod("show", "summary.cogarch.gmm",
          function (object)
          {
            
            cat("GMM estimation\n\nCall:\n")
            print(object@call)
            cat("\nCoefficients:\n")
            print(coef(object))
            cat("\n",paste0(paste("Log.objFun", object@objFun),":"), object@objFunVal, "\n")
            #cat("objFun", object@min, "\n")
          }
)


setMethod("summary", "cogarch.gmm.incr",
          function (object, ...)
          {
            cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
            m2logL <- 0
            data<-Re(coredata(object@Incr.Lev))
            data<- data[!is.na(data)]
            obj<- object@min
            labFun <- object@objFun
            
            tmp <- new("summary.cogarch.gmm.incr", call = object@call, coef = cmat,
                       m2logL = m2logL,
                       objFun = labFun,
                       objFunVal = obj,
                       MeanI = mean(data),
                       SdI = sd(data),
                       logLI = object@logL.Incr,
                       TypeI = object@model@measure.type,
                       NumbI = length(data),
                       StatI =summary(data)
            )
            tmp
          }
)
# 
setMethod("show", "summary.cogarch.gmm.incr",
          function (object)
          {
            
            cat("Two Stages GMM estimation \n\nCall:\n")
            print(object@call)
            cat("\nCoefficients:\n")
            print(coef(object))
#             cat("\n-2 log L:", object@m2logL, "\n")
            cat("\n",paste0(paste("Log.objFun", object@objFun),":"), object@objFunVal, "\n")
            
            cat(sprintf("\n\nNumber of increments: %d\n",object@NumbI))
            cat(sprintf("\nAverage of increments: %f\n",object@MeanI))
            cat(sprintf("\nStandard Dev. of increments: %f\n",object@SdI))
            if(!is.null(object@logLI)){
              
              cat(sprintf("\n\n-2 log L of increments: %f\n",-2*object@logLI))
            }
            cat("\nSummary statistics for increments:\n")
            print(object@StatI)
            cat("\n")
          }
)


