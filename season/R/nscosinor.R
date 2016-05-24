##**********************************
##** Seasonal decomposition macro **
##**********************************
##** Adrian Barnett               **
##** April 2008, updated Dec 2011 **
##**********************************
### Inputs
## data = data
## cycles = cycles (e.g. in months, f=c(6,12))
## tau = vector of smoothing parameters, tau[1] for trend, tau[2] for 1st seasonal parameter, tau[2] for 2nd seasonal parameter, etc
## niters = total number of MCMC samples (default=1000)
## burnin = number of MCMC samples discarded as a burn-in (default=500)
## lambda = distance between observations (lambda=1/12 for monthly data)
## div = divisor at which MCMC sample progress is reported (default=50)
## monthly = TRUE for monthly data
### Outputs
## trend = mean trend and 95% confidence interval 
## season = mean season and 95% confidence interval 
## residuals = residuals (based on mean trend and season)
## stats = estimated amplitude, phases and noise
## chains = MCMC chain of variance estimates (std.error for overall sd(error), std.season for seasonal parts)
##
## assumes year and month exist in data; assumes no missing data
`nscosinor` <-
  function(data,response,cycles,niters=1000,burnin=500,tau,
           lambda=1/12,div=50,monthly=TRUE,alpha=0.05){
    names<-names(data)
    yearyes<-sum(names=='year')
    monthyes<-sum(names=='month')
    if (yearyes<1|monthyes<1) {
      stop("Data needs to contain numeric year and month variables")}
    if (length(tau)!=length(cycles)+1) {
      stop("Need to give a smoothing parameter (tau) for each cycle, plus one for the trend")}
    resp=subset(data,select=response)[,1] # instead of attach
    if (sum(is.na(resp))>0) {
      stop("Missing data in the dependent variable not allowed")}
    if (sum(cycles<=0)>0) {stop("Cycles cannot be <=0")}
    if (burnin>niters) {
      stop("Number of iterations must be greater than burn-in")}
###  was    yrmon<-year+((month-1)/12)
    yrmon<-data$year+((data$month-1)/12)  # 
    n<-length(resp);
    k<-length(cycles);
    kk<-2*(k+1);
    ## Get initial values 
    good.inits=nscosinor.initial(data=data,response=response,lambda=lambda,tau=tau,n.season=k)
    vartheta<-sqrt(good.inits[1]) # Initial estimates of var theta
    w<-vector(length=k,mode="numeric")
    for (index in 1:k){ 
      w[index]<-good.inits[2] # Initial estimate of lambda (w)
    }
    ## Empty chain matrices and assign initial values
    ampchain<-matrix(0,niters+1,k)
    phasechain<-matrix(0,niters+1,k)
    alphachain<-array(0,c(kk,n+1,niters));
    varthetachain<-matrix(0,niters+1)
    lchain<-matrix(0,niters+1,k)
    lchain[1,]<-w
    varthetachain[1]<-vartheta
    cmean<-rep(10,kk) # starting value for C_j
    for (iter in 1:niters){ 
      result<-kalfil(resp,f=cycles,vartheta=varthetachain[iter], # changed response to resp
                     w=lchain[iter,],tau=tau,lambda=lambda,cmean=cmean)
      varthetachain[iter+1]<-result$vartheta 
      lchain[iter+1,]<-result$w 
      alphachain[,,iter]<-result$alpha 
      ampchain[iter+1,]<-result$amp
      phasechain[iter+1,]<-result$phase
      cmean<-result$cmean
      ## Output iteration progress 
      if (iter%%div==0){cat("Iteration number",iter,"of",niters,"\n",sep=" ")}
    }
    ## Get mean and percentiles of alpha (trend and season), and overall fitted values
    trend<-as.data.frame(matrix(0,n,3))
    season<-as.data.frame(matrix(0,n,3*k))
    oseason<-as.data.frame(matrix(0,n,3))
    new.fitted<-as.data.frame(matrix(0,n,3))
    names(trend)<-c('mean','lower','upper')
    names(oseason)<-c('mean','lower','upper')
    names(new.fitted)<-c('mean','lower','upper')
    allseasons<-matrix(data=NA,ncol=niters-burnin,nrow=n)
    snums<-((1:k)*2)+1
    for (i in 1:n){ 
      for (j in (burnin+1):niters){ 
        allseasons[i,j-burnin]<-sum(alphachain[snums,i,j])
      }
    }
    lprob=alpha/2;
    uprob=1-(alpha/2);
    lnum<-round((niters-burnin)*lprob);
    unum<-round((niters-burnin)*uprob);
    for.fitted=allseasons+alphachain[1, 1:n,(burnin+1):niters]
    for (i in 1:n){ 
      trend$mean[i]<-mean(alphachain[1,i,burnin:niters])
      trend$lower[i]<-sum(as.numeric(rank(alphachain[1,i,burnin:niters])==lnum)*alphachain[1,i,burnin:niters])
      trend$upper[i]<-sum(as.numeric(rank(alphachain[1,i,burnin:niters])==unum)*alphachain[1,i,burnin:niters])
      for (j in 2:(k+1)){ 
        snum<-((j-1)*2)+1
        season[i,((j-1)*3)-2]<-mean(alphachain[snum,i,burnin:niters])
        season[i,((j-1)*3)-1]<-sum(as.numeric(rank(alphachain[snum,i,burnin:niters])==lnum)*alphachain[snum,i,burnin:niters])
        season[i,((j-1)*3)]<-sum(as.numeric(rank(alphachain[snum,i,burnin:niters])==unum)*alphachain[snum,i,burnin:niters])
      }
      ## overall season
      oseason$mean[i]<-mean(allseasons[i,])
      oseason$lower[i]<-sum(as.numeric(rank(allseasons[i,])==lnum)*allseasons[i,])
      oseason$upper[i]<-sum(as.numeric(rank(allseasons[i,])==unum)*allseasons[i,])
      ## fitted values (with CIs)
      new.fitted$mean[i]<-mean(for.fitted[i])
      new.fitted$lower[i]<-sum(as.numeric(rank(for.fitted[i,])==lnum)*for.fitted[i,])
      new.fitted$upper[i]<-sum(as.numeric(rank(for.fitted[i,])==unum)*for.fitted[i,])
    }
    names(season)<-rep(c('mean','lower','upper'),k)
    ## Time
    if (monthly==TRUE){time<-yrmon}
    if (monthly!=TRUE){time<-1:n}
    ## Calculated fitted values and residuals
    fitted<-trend$mean+oseason$mean
    res<-resp-fitted # calculate the residuals
    ## original call with defaults (see amer package)
    ans <- as.list(match.call())
    frmls <- formals(deparse(ans[[1]]))
    add <- which(!(names(frmls) %in% names(ans)))
    call<-as.call(c(ans, frmls[add]))
    ## Returns
    toret<-list()
    toret$call<-call 
    toret$time<-time
    toret$trend<-trend
    toret$season<-season
    toret$oseason<-oseason
#    toret$fitted.values<-fitted # over-ride with values below
    toret$fitted.values<-new.fitted
    toret$residuals<-res
    toret$n<-n
    toret$chains$std.error<-varthetachain
    toret$chains$std.season<-matrix(data=NA,nrow=niters+1,ncol=k)
    toret$chains$phase<-matrix(data=NA,nrow=niters+1,ncol=k)
    toret$chains$amplitude<-matrix(data=NA,nrow=niters+1,ncol=k)
    for (i in 1:k){ # for multiple cycles
      toret$chains$std.season[,i]<-lchain[,i]
      toret$chains$phase[,i]<-phasechain[,i]
      toret$chains$amplitude[,i]<-ampchain[,i]
    }
    toret$chains <- coda::mcmc(cbind(
      toret$chains$std.error[(burnin+2):(niters + 1)], 
      toret$chains$std.season[(burnin+2):(niters +1), ],
      toret$chains$phase[(burnin+2):(niters + 1), ],
      toret$chains$amplitude[(burnin+2):(niters + 1), ]), start = burnin+1)
    # Add the names
    colnames(toret$chains)=c('std.error',paste('std.season',1:k,sep=''),paste('phase',1:k,sep=''),paste('amplitude',1:k,sep=''))
    toret$cycles<-cycles
    class(toret)<-'nsCosinor'
    return(toret)
  } # end of function

