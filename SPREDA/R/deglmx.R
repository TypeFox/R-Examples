deglmx <-
function(fixed, data, dyn.data, id, time, random=NULL, linear=TRUE, ytrend, splinetrend=NULL, splinetype=NULL,
                degree=NULL, knots=NULL, weights=NULL, subset=NULL, start, maxiter=10, method='BFGS', ... ){
  mfun <- match.call(expand.dots = FALSE)
  miss.data <- missing(data) || !is.data.frame(data)
  if (!inherits(fixed, "formula") || length(fixed) != 3) {
    stop("\nfixed-effects model must be a formula of the form \"resp ~ pred\"")
  }
  if(missing(dyn.data)){
    if(missing(random)){
      res=lm(fixed, data,  ...)
    }else{
      res=lme(fixed, data, random, weights, subset, ...)
    } 
    return(res)
  }else{
    y.name=deparse(fixed[[2]])
    covnames.tmp=getCov(fixed) # it can include constant covariate which will be added later in the item of data
    #extract fix covariates
    fixcov.index=match(covnames.tmp, colnames(data), 0L)
    fixcovnames=NULL
    if(sum(fixcov.index)>0){
      fixcovnames=colnames(data[, fixcov.index])
    }
    
    dyncov.index=match(covnames.tmp, colnames(dyn.data), 0L)
    dyncovnames=colnames(dyn.data[, dyncov.index])
    dyn.cov=dyn.data[, dyncov.index]  # dynamic covariate matrix
    rownames(dyn.cov)=NULL
    n.dyn=length(degree)
    
    check.index=length(unique(c(length(degree), length(knots), length(splinetrend), length(splinetype), length(dyncovnames))))
    if(check.index!=1) print("Please make sure the number of dynamic covariates are consistent")
    
    cov.mat=df=lam=NULL

   
    for(i in 1:n.dyn){
      tmp=splinetrend[i]*MIC.splines.basis(dyn.cov[,i], df = knots[i]+degree[i], knots = NULL,
                                  boundary.knots=NULL,type=splinetype[i], degree = degree[i], delta=0.01, eq.alloc=FALSE)$mat
      df.tmp=ncol(tmp)
      df=c(df, df.tmp)
      colnames(tmp)=paste(dyncovnames[i], 1:df.tmp, sep="")
      cov.mat=cbind(cov.mat, tmp)
      lam.tmp=rep(1, df.tmp)
      if(splinetype[i]=="Cs"){lam.tmp[1]=0}
      lam=c(lam, lam.tmp)
    }
    cov.mat.tmp=as.data.frame(cbind(dyn.data[,c(id, time)], cov.mat)) ###important data for dynamic covariate
    cov.cum=xmat.to.cumsum(dat=cov.mat.tmp)
    cov.cum=as.data.frame(cbind(cov.mat.tmp[,c(1,2)], cov.cum))
    colnames(cov.cum)[c(1,2)]=colnames(cov.mat.tmp)[c(1,2)]
    

    newformula=formula(paste(y.name, "~", paste(c(id, time, fixcovnames),  collapse="+"), sep=""))
    
        
    mfix <- match(c("fixed", "data", "weights", "subset"), names(mfun), 0L)
    if (mfix[1]==0) stop("a formula argument is required")
    mfix.fun <- mfun[c(1, mfix)]  
    mfix.fun[[1]] <- as.name("model.frame")
    names(mfix.fun)[2] <- "formula"
    mfix.fun[[2]]=newformula 
    fixdat <- eval(mfix.fun, parent.frame())
    n.fixdat=ncol(fixdat)
    
    ### data matrix
    newdata.tmp=xmat.obj.to.xmat(fixdat, cov.cum, id, time)
    newdata=as.data.frame(cbind(newdata.tmp[, id], newdata.tmp[, time], newdata.tmp[,y.name], 
                                INT=1, fixdat[,-c(1,2,3)], time.b=ytrend*newdata.tmp[, time],
                                newdata.tmp[, -c(1:n.fixdat)])) 
    colnames(newdata)=c(id, time, y.name, "INT", fixcovnames, "time.b", colnames(newdata.tmp)[-c(1:n.fixdat)])
    cov.mat.tmp=match.dat.fun(newdata, cov.mat.tmp)
    dyn.data=match.dat.fun(newdata, dyn.data)
    ### random part
    random.names=NULL
    if(linear){
      
      if(length(random)!=0){
        rand.tmp=getRanName(random)
        if(rand.tmp$rhs!=id){print("Error:The conditional part in random should be equal to id")}
        int.index=match("-1", rand.tmp$lhs, 0L)
        if(int.index>0){random.names=rand.tmp$lhs[-int.index]}else{
          int.index2=match("1", rand.tmp$lhs, 0L)
          if(int.index2>0){random.names=c("INT", rand.tmp$lhs[-int.index2])}else{
            random.names=c("INT", rand.tmp$lhs)
          }
        }   
      }
      lam=c(0, rep(0,length(fixcovnames)), 1, lam) #add int and time.b
      dat.obj=list(dat=newdata, lam=lam, dfs=df)
      beta.index=as.logical(c(0, rep(0,length(fixcovnames)), 1, rep(1, sum(df)), rep(0, 2*length(random.names)))) #extract beta for dyn.cov: time.b+dyn.cov
      fit=clme(dat.obj=dat.obj,theta0=start,random.names=random.names, no.iter=maxiter, trace=F, method=method)

    }else{
    
      newdata=newdata[, -4] # id time y time.b dyn.cov
      newdata[,4]=-newdata[,4]
      lam=c(1, lam) #add time.b
      dat.obj=list(dat=newdata, lam=lam, dfs=df)
      coef=coefinitial.ftn(dat=dat.obj, start=start)
      n.ids=length(unique(newdata[,1]))
      random.eff=matrix(0, nrow=n.ids, ncol=2)
      beta.index=as.logical(c(0,0, rep(1, length(lam)), rep(0, 4)))
      fit=optim.ftn.2(dat=dat.obj, coef=coef, random.eff=random.eff,  iter=maxiter)
    }
    if(linear){type="linear"}else{
      type="nonlinear"
    }
    res=list(type=type, fit=fit, dat=dat.obj, dyn.mat=cov.mat.tmp, ytrend=ytrend, 
             dyncovnames=dyncovnames,dyn.data=dyn.data, beta.index=beta.index, call=mfun)
    
    class(res)="deglmx"
    return(res)
    
  }
  
}
