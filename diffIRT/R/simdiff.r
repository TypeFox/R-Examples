simdiff=function(N,nit,ai=NULL,vi=NULL,gamma=NULL,theta=NULL,ter=NULL,model="D",max.iter=19999,eps=1e-15){
  vp=theta
  ap=gamma
   if (missing(N)|missing(nit)) stop("Number of subjects and/or the number of items is missing.")
   
  if(toupper(model)!="D" & toupper(model)!="Q" )
    stop("The 'model' argument should be either 'D' for the D-diffusion IRT model or 'Q' for the Q-diffusion IRT model.")

  if(!is.null(ai) & !is.null(vi) & !is.null(ap) & !is.null(vp) & !is.null(ter)) {
    if(length(ai)!=length(vi) | length(ai)!=length(ter) | length(vi)!=length(ter) )
      stop("Error: 'a[i]', 'v[i]', and 'Ter[i]' are not the same length\n")
    if(length(ap)!=length(vp)) stop("Error: 'gamma[p]' and 'theta[p]' are not the same length\n")
    if(length(ai)!=nit | length(ap)!=N) 
      stop("Arguments 'N' and 'nit' are not in line with the length of the vectors of true parameter values.\n")
    if(sum((ter<0)*1)!=0)
      stop("Error: At least one of the true values for 'Ter[i]' is non-possitive.\n")
    if(sum((ai<0)*1)!=0)
      cat("Warning: true values for 'a[i]' are not strictly possitive which is not in line with the D- or Q-diffusion IRT model\n")
    if(sum((ap<0)*1)!=0)
      cat("Warning: true values for 'gamma[p]' are not strictly possitive which is not in line with the D- or Q-diffusion IRT model\n")
    if(model=="Q" & sum((vi<0)*1)!=0)
      cat("Warning: true values for 'v[i]' are not strictly possitive which is not in line with the Q-diffusion IRT model\n")
    if(model=="Q" & sum((vp<0)*1)!=0)
      cat("Warning: true values for 'theta[p]' are not strictly possitive which is not in line with the Q-diffusion IRT model\n")
  }
  else{
    if(is.null(ai) & is.null(vi) & is.null(ap) & is.null(vp) & is.null(ter)){
      ap=rlnorm(N,0,.3)
      ai=exp(runif(nit,-1,0))
      ter=runif(nit,.5,1.5)
      if(toupper(model)=="Q"){
        vp=rlnorm(N,0,.3)
        vi=exp(runif(nit,0,1))
      }
      if(toupper(model)=="D"){
        vp=rnorm(N,0,.5)
        vi=runif(nit,-1.5,.5)
      }
    } else stop("Error: true values for either a[i], v[i], gamma[p], theta[p], and/or Ter[p] are missing")
  }

  rt=x=matrix(,N,nit)

  for(j in 1:nit)
  for(p in 1:N){{
    a=ap[p]/ai[j]
    if(toupper(model)=="D") mu=vp[p]-vi[j]
    else mu=vp[p]/vi[j]
    si=1
    M=pi*si^2/a^2 * (exp(a*mu/(2*si^2))+exp(-a*mu/(2*si^2))) * 1/ (mu^2/(2*si^2)+pi^2*si^2 / (2*a^2))
    lmb = mu^2/(2*si^2) +pi^2*si^2/(2*a^2)
    ou=c(-999)
    rej=0
    numb=1
    while(ou==-999){
      v=runif(1); u=runif(1);
      FF=pi^2*si^4 * 1/(pi^2*si^4+mu^2*a^2)
      sh1=1
      sh2=0
      sh3=0
      i=0
      while(abs(sh1-sh2)>eps | abs(sh2-sh3)>eps){
        sh1=sh2
        sh2=sh3
        i=i+1
        sh3= sh2 + (2*i+1)*(-1)^i*(1-u)^(FF*(2*i+1)^2)
      }
    eval=1+(1-u)^-FF * sh3
    if(v<=eval) ou=1/lmb*abs(log(1-u))
    else rej=rej+1
    if(rej==max.iter) stop("Rejection algorithm failed. Increase the 'max.iter' argument or try different true parameter values.")
    }
    score=((1/(1+exp(-mu*a))) > runif(1))*1
    rt[p,j]=ou+ter[j]
    x[p,j]=score
  }
  }
  return(list(rt=rt,x=x,ai=ai,vi=vi,gamma=ap,theta=vp,ter=ter))
}

