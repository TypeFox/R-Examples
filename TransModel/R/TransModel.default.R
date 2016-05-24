TransModel.default <-
function(formula = formula(data),data = parent.frame(),r=0,CICB.st=FALSE,dx=0.001,iter.max=100,num.sim=200){
    copy.data <- data
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "copy.data", "na.action")
    mf <- mf[match(temp, names(mf), nomatch = 0)]
    mf[[1]] <- as.name("model.frame")
    temp.x<-terms(formula, data = copy.data)
    mf$formula <- temp.x
    mf <- eval(mf, parent.frame())
    Y <- model.extract(mf, "response")
    obs.t<-Y[,1]
    delta<-Y[,2]
    attr(temp.x, "intercept") <- 0
    
    z <- model.matrix(temp.x,copy.data)
    zc = t(t(z) - apply(z,2,mean))
    ix = order(obs.t)
    ord.delta = delta[ix]
    ord.z = as.matrix(zc[ix,])
    colnames(ord.z)= colnames(z)
    p = dim(z)[2]
    n = dim(z)[1]

    beta.ini = rep(0,p)
    Rt.ini = cumsum(1/(n:1))

    ###solving estimating equations for beta and Rt
    wt<-rep(1,n)
    temp.b<-solve.beta(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max)
    beta<-temp.b$Beta
    Rt<-temp.b$data$Rt
    converged<-temp.b$converged
    iter<-temp.b$iter
    ###compute the variance-covariance matrix of beta
    dRt = Rt - c(0,Rt[1:(n-1)])
    ord.bz = as.numeric(ord.z%*%beta)
    a = Rt*exp(-ord.bz)/(1+r*Rt*exp(-ord.bz))
    B2 = sum(a) - cumsum(a) + a
    b = exp(-ord.bz)/(1+r*Rt*exp(-ord.bz))^2
    B1 = sum(b) - cumsum(b) + b
    bt = exp(cumsum(dRt*B1/B2))
    Bzt = matrix(0,nrow=n,ncol=p)
    Zbart = matrix(0,nrow=n,ncol=p)
    A = matrix(0,nrow=p,ncol=p)
    V = matrix(0,nrow=p,ncol=p)
    for(i in 1:(n-1)){
      B = rep(bt[i],n-i+1)/bt[i:n]
      if(p>1) Bzt[i,] = apply(ord.z[i:n,]*a[i:n]*B,2,sum)
	if(p==1) Bzt[i,] =sum(ord.z[i:n,]*a[i:n]*B)
      Zbart[i,] = Bzt[i,]/B2[i]
      c1 = exp(-ord.bz)*dRt[i]/(1+r*Rt[i]*exp(-ord.bz))^2
      A = A + (t(ord.z[i:n,]) - as.numeric(Zbart[i,]))%*%(ord.z[i:n,]*c1[i:n])
      c2 = exp(-ord.bz)*dRt[i]/(1+r*Rt[i]*exp(-ord.bz))
      V = V + (t(ord.z[i:n,]) - as.numeric(Zbart[i,]))%*%(t(t(ord.z[i:n,]) - as.numeric(Zbart[i,]))*c2[i:n])
    }
    Bzt[n,] = ord.z[n,]*a[n]
    Zbart[n,] = Bzt[n,]/B2[n]
    d = as.numeric(ord.z[n,] - Zbart[n,])
    A = A + d%*%t(as.numeric(ord.z[n,])*exp(-ord.bz[n])*dRt[n]/(1+r*Rt[n]*exp(-ord.bz[n]))^2)
    V = V + d%*%t(d*exp(-ord.bz[n])*dRt[n]/(1+r*Rt[n]*exp(-ord.bz[n])))
    invA = solve(A)
    Sigma = invA%*%V%*%t(invA)
           
    vcov=matrix(Sigma,nrow=p)
    names(beta)<-colnames(vcov)<-rownames(vcov)<-all.vars(formula)[-c(1:2)]

    output<-list(coefficients=-beta,vcov=vcov,formula=formula,nz=ncol(z),CICB.st=CICB.st,
         r=r,converged=c(converged,iter),df=n-1,call=match.call())

    if(CICB.st){
      cicb<-var.e(obs.t,delta,z,r,dx,iter.max,n,p,num.sim)
      output$rep.beta<-cicb[[1]]
      output$rep.Rt<-cicb[[2]]
    }

    output$ord.time<-sort(obs.t)
    output$delta=ord.delta
    output$z=z
    output$Rt=temp.b$data$Rt
    class(output)<-"TransModel"
    return(output)
}