GQD.estimates = function(x,thin = 100, burns, CI = c(0.05,0.95), corrmat = FALSE, acf.plot =TRUE)
{
  if(class(x)=='GQD.mle')
  {
    sigma = sqrt(diag(solve(-x$opt$hessian)))
    upper = x$opt$par+1.96*sigma
    lower = x$opt$par-1.96*sigma
    EstCI = data.frame(Estimate=x$opt$par, Lower_95=lower,Upper_95=upper)
    form  = function(x,mm = 2){format(round(x, mm), nsmall = mm)}
    rownames(EstCI) <- paste0('theta[',1:length(x$opt$par),']')
    dat2 =data.frame(form(solve(-x$opt$hessian)/(sigma%o%sigma),2))
    rownames(dat2) <- paste0('theta[',1:length(x$opt$par),']')
    colnames(dat2) <- paste0('theta[',1:length(x$opt$par),']')
    if(corrmat){return(list(estimates=data.frame(form(EstCI,3)),corrmat = dat2))}
    return(data.frame(round(EstCI,3)))
  }
  
  if(class(x)=='GQD.mcmc')
  {
    if(missing(burns)){burns =min(round(dim(x$par.matrix)[1]/2),25000)}
    windw = seq(burns,dim(x$par.matrix)[1],thin)
    est = apply(x$par.matrix[windw,], 2, mean)
    CI=t(apply(x$par.matrix[windw,], 2, quantile,probs = CI))
    form = function(x,mm = 2){format(round(x, mm), nsmall = mm)}
    dat=data.frame(cbind(form(cbind(est,CI),3)))
    rownames(dat)=paste0('theta[',1:dim(x$par.matrix)[2],']')
    colnames(dat) = c('Estimate','Lower_CI','Upper_CI')
    
    dat2=data.frame(form(cor(x$par.matrix[windw,])))
    rownames(dat2)=paste0('theta[',1:dim(x$par.matrix)[2],']')
    colnames(dat2)=paste0('theta[',1:dim(x$par.matrix)[2],']')
    if(acf.plot)
    {
      nper=dim(x$par.matrix)[2]
      if(nper==1){par(mfrow=c(1,2))}
      if(nper==2){par(mfrow=c(2,2))}
      if(nper==3){par(mfrow=c(2,2))}
      if(nper>3)
      {
        d1=1:((nper)+1)
        d2=d1
        O=outer(d1,d2)
        test=O-((nper)+1)
        test[test<0]=100
        test=test[1:4,1:4]
        test
        wh=which(test==min(test))
        
        d1=d1[col(test)[wh[1]]]
        d2=d2[row(test)[wh[1]]]
        par(mfrow=c(d1,d2))
      }
      cols=rainbow_hcl(nper, start = 10, end = 275,c=100,l=70)
      for(i in 1:dim(x$par.matrix)[2])
      {
        acf(x$par.matrix[windw,i],main=paste0('ACF: theta[',i,']\nThin=',thin,', Burns=',burns,', N=',length(windw)),col = cols[i],lwd=2)
      }
    }
    if(corrmat){return(list(estimates = dat, corrmat = dat2))}
    return(dat)
  }
  
}