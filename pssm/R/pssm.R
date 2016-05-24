pssm <-
  function(progr,survv,dat,intervals=5,start=NULL,rescale=1){
    #avoid note on lack of declaration for variables declared by "assign"
    loglike=rescale1=m=m1=outp=message=NULL
  prior=NULL #prior is not used in estimation but retained in log-likelihood
  #for metropolis hastings
	call=match.call()	
    #library(MASS)
    m=intervals
    tr1<-function(x) tryCatch(x,error=function(e) stop('issue with input data frame',call.=FALSE))
    tr2<-function(x) tryCatch(x,error=function(e) {return(NULL)})
    # Function to perform optimization
	sp=is.null(progr)
    ss=is.null(survv)
    both=!(ss|sp)
    cov1=c()
    cov2=c()
	#extracts covariates
    fr=function(a,b) if(b<a) NULL else a:b  #Need to fix this wherever it is used
  error=NULL
  namesdat=names(dat)
  tst<-function(x){!all(is.element(x,namesdat))}
  tst1<-function(nx){
    x=dat[,nx]
    mx=max(x);
    mn=min(x);
    if(is.na(mx)) return(TRUE) 
    else {if (mn<0) return(TRUE) 
    else return(FALSE)}
  }
  rtn<-function(){print("Database doesn't have all the call variables");
                  assign("error","ERROR",parent.frame())
                  return(NULL)}
  rtn1<-function(){print("Problem with the call variables");
                   assign("error","ERROR",parent.frame())
                   return(NULL)}
  if(!sp){
      prog<-as.list(attr(terms(progr),"variables"))
      tp0<-as.character(prog[2][[1]][[2]])
      tp1<-as.character(prog[2][[1]][[3]])
      cov1<-as.character(prog[-(1:2)])
      if (tst(c(tp0,tp1,cov1))) {rtn()}
      if (tst1(tp0)) rtn1()
      if (any(ifelse(!is.na(dat[,tp1]),dat[,tp1]<dat[,tp0],FALSE))) rtn1()
      if (any(c(is.na(dat[,cov1])))) rtn1()
      #error=function(x) print("Problem in call for progression model"))
    }
    
    if(!ss){
      #Progression in Model
      surv=as.list(attr(terms(survv),"variables"))
      cov2=as.character(surv[-(1:2)])
      td=as.character(surv[2][[1]][[2]])
      cp=as.character(surv[2][[1]][[3]])
      if (tst(c(td,cp,cov2))) {rtn()}
      if (tst1(td)) rtn1()
      if (tst1(cp)) rtn1()
      if (any(c(is.na(dat[,cov1])))) rtn1()
    }
   
    if(!is.null(error)) {print(error);return(NULL)}
	   if (!(sp||ss)) {tmax=max(c(dat[,tp1],dat[,td]),na.rm=TRUE)}
	   else {if (!ss) tmax=max(dat[,td],na.rm=TRUE) else tmax=max(dat[,tp1],na.rm=TRUE)}  
     if(tmax>m) rescale=m/tmax
  
  opt<-function(funct){
		m1=m
        rescale1=rescale
        message=""
		fname=match.call()[[2]]
		mcovs=length(cov1)+length(cov2)
		if (is.null(prior)) prior=rep(0,mcovs)
    restrt=function(x){

    ms1=m1+m1*(m1+1)/2+mcovs
    ms2=m1+mcovs
    mc=length(x)
    if ((mc==ms1)||(mc==ms2)) return(x) else {
      if(mc==(ms2+1)) return(c(rep(mean(x[1:(m1+1)]),m1),x[(m1+1):(m1+mcovs)])) else 
      return(c(rep(mean(x[1:((m1+2)*(m1+1)/2)]),(m1*(m1+1)/2)),
               rep(mean(x[(((m1+2)*(m1+1)/2)+1):(((m1+2)*(m1+1)/2)+m1)]),m1),
               x[(((m1+2)*(m1+1)/2)+m1+1):(((m1+2)*(m1+1)/2)+m1+mcovs)]))}
		}
    
      while(TRUE){
         if(fname=="llikef"){
		    if(is.null(start)) st=c(startv(dt,m1),rep(0,length(cov1)+length(cov2))) else st=restrt(start)
	        fc<-funct(cov1,cov2,dt,m=m1,accumulate=TRUE,rescale=rescale1,prior=prior)}
		else {
			if(fname=="rprog"){
			if(is.null(start)) st=c(startv(dt,intervals,method='progression'),rep(0,length(cov1))) else st=restrt(start)
			   fc<-funct(cov1,dt,m=m1,accumulate=TRUE,rescale=rescale1)}
			else {
			if(is.null(start)) st=c(startv(dt,m1,method='survival'),rep(0,length(cov2))) else st=restrt(start)
			  fc<-funct(cov2,dt,m=m1,accumulate=TRUE,rescale=rescale1)}
			}	
        outp=tr2(optim(st,fc, method="BFGS",hessian=TRUE,control=list(fnscale=-1)))
        if(is.null(outp)||(outp$convergence>0)) outp=tr2(optim(st,fc, method="Nelder-Mead",hessian=TRUE,control=list(fnscale=-1,maxit=1000)))  else break
        if(outp$convergence>0){
          if(m1>1){
            message="Number of Intervals Reduced"
            m1=m1-1
            rescale1=rescale1*m1/(m1+1)} else {
              message="optimization failed"
              break        
            }
        } else break
      }
      assign("loglike",fc,envir=parent.frame())
      assign("rescale1",rescale1,envir=parent.frame())
      assign("m",m1,envir=parent.frame())
      assign("outp",outp,envir=parent.frame())
      assign("message",message,envir=parent.frame())
      return(outp$convergence)
    }
	 if(!(ss|sp)){
      #Full Model
      dt=tr1(data.frame(tprog0=dat[,tp0],tprog1=dat[,tp1],tdeath=dat[,td],cdeath=dat[,cp],dat[,cov1],dat[,cov2]))
  
  	  
	  names(dt)<-c("tprog0","tprog1","tdeath","cdeath",cov1,cov2)


      convergence=opt(llikef)
  	 
      mms=m+m*(m+1)/2
      ep=fr((mms+1),(mms+length(cov1)))
      ehp=fr((m*(m+1)/2+1),mms)
      es=fr((mms+length(cov1)+1),(mms+length(cov1)+length(cov2)))
      ehs=fr(1,(m*(m+1)/2))
    } else {
      if(!ss){
        #Has only Survival
        mms=m   
        dt=tr1(data.frame(tdeath=dat[,td],cdeath=dat[,cp],dat[,cov2]))

        names(dt)<-c("tdeath","cdeath",cov2)
        
        convergence=opt(rsurv)
        ep=quote(NULL)
        ehp=quote(NULL)
        es=fr(mms+length(cov1)+1,mms+length(cov1)+length(cov2))
        ehs=fr(1,m)}
      else {
        #Has only Progression
        mms=m
        ep=fr((mms+1),(mms+length(cov1)))
        ehp=quote(1:m)
        es=quote(NULL)
        ehs=quote(NULL)
        es=quote(NULL)
        ehs=quote(NULL)
        dt=tr1(data.frame(tprog0=dat[,tp0],tprog1=dat[,tp1],dat[,cov1]))
	 
        names(dt)<-c("tprog0","tprog1",cov1)
        
        convergence=opt(rprog)
        mms=m
        ep=fr((mms+1),(mms+length(cov1)))
        ehp=quote(1:m)
        es=quote(NULL)
        ehs=quote(NULL)
        es=quote(NULL)
        ehs=quote(NULL)
      }	 
    }	
    if(!is.null(outp))
    {
      cov=ginv(-outp$hessian)
      dcov=diag(cov)
      if(any(dcov<=0)) print("Hessian is not positive definite, consider reducing the number of intervals")
      se=sqrt(pmax(diag(cov),0))
      
      outpssm=new("pssm",
		          call=call,
                  convergence=0,
                  loglike=loglike,
                  estimates=outp$par,
                  se.estimates=se,
                  covariance.estimates=cov,
                  estimates.progression=outp$par[eval(ep)],
                  se.estimates.progression=se[eval(ep)],
                  estimates.survival=outp$par[eval(es)],
                  se.estimates.survival=se[eval(es)],
                  hazard.progression=outp$par[eval(ehp)],
                  hazard.survival=outp$par[eval(ehs)],
                  intervals=as.integer(m),
                  rescale=rescale1,
                  formula.progression=as.formula(progr),
                  formula.survival=as.formula(survv),
                  progression.covariate.list=as.character(cov1),
                  survival.covariate.list=as.character(cov2),message=message)
    }
    else outpssm=NULL
    return(outpssm)
  }
