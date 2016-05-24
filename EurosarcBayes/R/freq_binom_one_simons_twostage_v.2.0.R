


freq_binom_one_simons_twostage=function(p0,p1,alpha,power,prior.a=0,prior.b=0,nmax=100,round=TRUE,method="optimal"){
  # obtain solutions from ph2simon
  solutions=ph2simon(p0,p1,alpha,1-power,nmax)

  # find solutions
  opt=which.min(solutions$out[,5])
  opt=solutions$out[opt,]
  minimax=solutions$out[1,]


  if(method=="optimal"){
    return(properties_binom_one(failure=opt[c(1,3)],success=as.numeric(c(NA,opt[3]+1)),reviews=opt[c(2,4)],p0,p1,prior.a,prior.b,round))
  } else if(method=="minmax"){
    return(properties_binom_one(failure=minimax[c(1,3)],success=c(NA,minimax[3]+1),reviews=minimax[c(2,4)],p0,p1,prior.a,prior.b,round))
  }

}

# ended
