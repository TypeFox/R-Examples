##### covariance matrix estimator

para=function(ivmodel){
  output=list()
  fit1=lm(ivmodel$Dadj~ivmodel$Zadj-1)
  output$gamma=coef(fit1)
  names(output$gamma)=NULL
  output$beta=c(KClass(ivmodel, k=1)$point.est)
  hat_eps=ivmodel$Yadj-c(output$beta)*ivmodel$Dadj
  hat_eta=resid(fit1)
  output$sigmau=sqrt(sum(hat_eps^2)/(ivmodel$n-ivmodel$p))
  output$sigmav=sqrt(sum(hat_eta^2)/(ivmodel$n-ivmodel$p))
  output$rho=sum(hat_eta*hat_eps)/(ivmodel$n-ivmodel$p)/output$sigmau/output$sigmav
  return(output)
}

#####  AR test and CI

AR.test=function(ivmodel, beta0=0, alpha=0.05){
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; 
  n=ivmodel$n;  k=ivmodel$p;  l=ivmodel$L

  if(ncol(Yadj)>1){
    print("The outcome variable should be in one dimension!")
	return(NULL)
  }
  if(ncol(Dadj)>1){
    print("The treatment variable should be in one dimension!")
	return(NULL)
  }
  if(l+k>=n){
    print("Too many IVs, AR can't handle!")
    return(NULL)
  }

### test  
  temp=Yadj-beta0*Dadj
  Fstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(Fstat, df1=l, df2=n-k-l)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l)*l/(n-k-l)
  coef.beta0sq=cval*sum(Dadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Dadj)^2)
  coef.beta0=-2*cval*sum(Dadj*Yadj)+2*(cval+1)*sum(Dadj*qr.fitted(ivmodel$ZadjQR, Yadj))
  coef.constant=cval*sum(Yadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Yadj)^2)
  Delta=coef.beta0^2-4*coef.constant*coef.beta0sq

  ci=matrix(NA, ncol=2)
  colnames(ci)<-c("lower", "upper")

  if(coef.beta0sq==0){
    if(coef.beta0>0){
      info=c("[",-coef.constant/coef.beta0,",Infinity)")
      ci[1,]=c(-coef.constant/coef.beta0, Inf)
    }
    if(coef.beta0<0){
      info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
      ci[1,]=c(-Inf, -coef.constant/coef.beta0)
    }
    if(coef.beta0==0){
      if(coef.constant>=0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.constant<0){
        info="Empty Set"
      }
    }
  }
  
  if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",lower.root,",",upper.root,"]")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }  

  return(list(Fstat=Fstat, df=c(l, n-k-l), p.value=p.value,
              ci.info=info, ci=ci))
}


AR.power=function(n, k, l, beta, gamma, Zadj_sq, 
                  sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
	print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*n*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	    (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  temp = qf(1-alpha, df1=l, df2=n-k-l)
  power = 1-pf(temp, df1=l, df2=n-k-l, ncp=ncp)

  return(power)
}

AR.size=function(power, k, l, beta, gamma, Zadj_sq, 
                 sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
    print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	  (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  oldn<-k+l+1
  state<-1
  while(state){
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*oldn)  
    if(temppower < power){
	  oldn <- oldn*2
	}else{
	  state <- 0
	}
  }

  lower <- oldn%/%2
  upper <- oldn
  while((upper-lower)>2){
    new <- (upper+lower)%/%2
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*new)  
    if(temppower < power){
	  lower <- new
	}else{
	  upper <- new
	}
  }

  return(upper)
}


#####  AR sensitivity test and CI

ARsens.test=function(ivmodel, beta0=0, alpha=0.05, deltarange=NULL){
  if(is.null(deltarange))
    return(NULL)

  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; 
  n=ivmodel$n;  k=ivmodel$p;  l=ivmodel$L

  if(ncol(Yadj)>1){
    print("The outcome variable should be in one dimension!")
	return(NULL)
  }
  if(ncol(Dadj)>1){
    print("The treatment variable should be in one dimension!")
	return(NULL)
  }
  if(l+k>=n){
    print("Too many IVs, AR can't handle!")
    return(NULL)
  }
  
  if(l!=1){
    print("Please input exact one IV for AR sensitivity analysis")
	return(NULL)
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    ncp=max(deltarange^2)*sum(Zadj^2)
  }else{
    print("Wrong input of the sensitivity range.")
    return(NULL)
  }
  
### test  
  temp=Yadj-beta0*Dadj
  ncFstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(ncFstat, df1=l, df2=n-k-l, ncp=ncp)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l, ncp=ncp)*l/(n-k-l)
  coef.beta0sq=cval*sum(Dadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Dadj)^2)
  coef.beta0=-2*cval*sum(Dadj*Yadj)+2*(cval+1)*sum(Dadj*qr.fitted(ivmodel$ZadjQR, Yadj))
  coef.constant=cval*sum(Yadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Yadj)^2)
  Delta=coef.beta0^2-4*coef.constant*coef.beta0sq

  ci=matrix(NA, ncol=2)
  colnames(ci)<-c("lower", "upper")

  if(coef.beta0sq==0){
    if(coef.beta0>0){
      info=c("[",-coef.constant/coef.beta0,",Infinity)")
      ci[1,]=c(-coef.constant/coef.beta0, Inf)
    }
    if(coef.beta0<0){
      info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
      ci[1,]=c(-Inf, -coef.constant/coef.beta0)
    }
    if(coef.beta0==0){
      if(coef.constant>=0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.constant<0){
        info="Empty Set"
      }
    }
  }
  
  if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",lower.root,",",upper.root,"]")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }  

  return(list(ncFstat=ncFstat, df=c(l, n-k-l), ncp=ncp, 
              p.value=p.value, ci.info=info, ci=ci, deltarange=deltarange))
}

ARsens.size=function(power, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                     rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  if(ncp1<=ncp2){
    print("Sensitivity range too large")
	stop()
  }else{
    oldn<-k+2
    state<-1
    while(state){
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*oldn)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*oldn)  
      if(temppower < power){
	    oldn <- oldn*2
	  }else{
	    state <- 0
	  }
    } 

    lower <- oldn%/%2
    upper <- oldn
    while((upper-lower)>2){
      new <- (upper+lower)%/%2
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*new)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*new)  
      if(temppower < power){
	    lower <- new
	  }else{
	    upper <- new
	  }
    }
  }
  return(upper)
}

#####  calculate the power of AR sensitivity analysis

ARsens.power=function(n, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                      rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*n*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     n*Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*n*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  temp = qf(1-alpha, df1=1, df2=n-k-1, ncp=ncp2)
  power = 1-pf(temp, df1=1, df2=n-k-1, ncp=ncp1)

  return(power)
}

TSLS.power=function(n, beta, rho_ZD, sigmau, sigmaDsq, alpha=0.05){
  return(c(1+pnorm(-qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)-pnorm(qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)))
}

TSLS.size=function(power, beta, rho_ZD, sigmau, sigmaDsq, alpha=0.05){
  return(c((qnorm(1-alpha/2)+qnorm(power))^2*sigmau^2/beta^2/rho_ZD^2/sigmaDsq))
}


# ### calculating power for TSLS AR and ARsens

# IV.power=function(n, alpha=0.05, beta, type="TSLS", sigmau,
                  # rho_ZD=NULL, sigmaDsq=NULL,
				  # k=NULL, l=NULL, gamma=NULL, Zadj_sq=NULL, sigmav=NULL, rho=NULL,
				  # deltarange=NULL, delta=NULL){
  # if(type=="TSLS" & !is.null(rho_ZD) & !is.null(sigmaDsq)){
    # return(c(1+pnorm(-qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)-pnorm(qnorm(1-alpha/2)-beta*rho_ZD*sqrt(n*sigmaDsq)/sigmau)))
  # }
  # if(type=="AR" & !is.null(k) & !is.null(l) & !is.null(gamma) & !is.null(Zadj_sq) & !is.null(sigmav) & !is.null(rho)){
    # if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
	  # print("The dimension of Zadj_sq doesn't match gamma")
      # stop()
    # }
    # ncp=beta^2*n*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	    # (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

    # temp = qf(1-alpha, df1=l, df2=n-k-l)
    # power = 1-pf(temp, df1=l, df2=n-k-l, ncp=ncp)

    # return(power)  
  # }

  # if(type=="ARsens" & !is.null(k) & !is.null(gamma) & !is.null(Zadj_sq) & !is.null(sigmav) & !is.null(rho) & !is.null(deltarange)){  
    # if(!is.numeric(gamma) | length(gamma)!=1){
	  # print("Wrong input of gamma or gamma is not one dimension")
	  # stop()
    # }
  
    # if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	  # print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	  # stop()
    # }

    # if(is.numeric(deltarange) & length(deltarange)==2){
      # deltarange<-sort(deltarange)
      # ncp2=max(deltarange^2)*n*Zadj_sq
	  # if(is.null(delta)){
	    # if((deltarange[1] < -gamma*beta/sigmau) & 
	       # (deltarange[2] > -gamma*beta/sigmau)){
	      # ncp1=0
	    # }else{
	      # ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         # (beta*gamma+deltarange[2]*sigmau)^2)*
		     # n*Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	    # }	
	  # }else{
	    # ncp1=(beta*gamma+delta*sigmau)^2*n*Zadj_sq/
	       # (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  # }
    # }else{
      # print("Wrong input of the sensitivity range.")
      # stop()
    # }

    # temp = qf(1-alpha, df1=1, df2=n-k-1, ncp=ncp2)
    # power = 1-pf(temp, df1=1, df2=n-k-1, ncp=ncp1)

    # return(power)  
  # }
  # print("Input error.")
  # return(NULL)
# }


# ### calculating size for TSLS AR and ARsens

# IV.size=function(power, alpha=0.05, beta, type="TSLS", sigmau,
                  # rho_ZD=NULL, sigmaDsq=NULL,
				  # k=NULL, l=NULL, gamma=NULL, Zadj_sq=NULL, sigmav=NULL, rho=NULL,
				  # deltarange=NULL, delta=NULL){
  # if(type=="TSLS" & !is.null(rho_ZD) & !is.null(sigmaDsq)){
    # return(c((qnorm(1-alpha/2)+qnorm(power))^2*sigmau^2/beta^2/rho_ZD^2/sigmaDsq))
  # }
  # if(type=="AR" & !is.null(k) & !is.null(l) & !is.null(gamma) & !is.null(Zadj_sq) & !is.null(sigmav) & !is.null(rho)){
    # if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
      # print("The dimension of Zadj_sq doesn't match gamma")
	  # stop()
    # }
    # ncp=beta^2*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	    # (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

    # oldn<-k+l+1
    # state<-1
    # while(state){
      # temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
      # temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*oldn)  
      # if(temppower < power){
	    # oldn <- oldn*2
	  # }else{
	    # state <- 0
	  # }
    # }

    # lower <- oldn%/%2
    # upper <- oldn
    # while((upper-lower)>2){
      # new <- (upper+lower)%/%2
      # temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
      # temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*new)  
      # if(temppower < power){
	    # lower <- new
	  # }else{
	    # upper <- new
	  # }
    # }

    # return(upper)    
  # }

  # if(type=="ARsens" & !is.null(k) & !is.null(gamma) & !is.null(Zadj_sq) & !is.null(sigmav) & !is.null(rho) & !is.null(deltarange)){  
    # if(!is.numeric(gamma) | length(gamma)!=1){
	  # print("Wrong input of gamma or gamma is not one dimension")
	  # stop()
    # }
  
    # if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	  # print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	  # stop()
    # }

    # if(is.numeric(deltarange) & length(deltarange)==2){
      # deltarange<-sort(deltarange)
      # ncp2=max(deltarange^2)*Zadj_sq
	  # if(is.null(delta)){
	    # if((deltarange[1] < -gamma*beta/sigmau) & 
	       # (deltarange[2] > -gamma*beta/sigmau)){
	      # ncp1=0
	    # }else{
	      # ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         # (beta*gamma+deltarange[2]*sigmau)^2)*
		     # Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	    # }	
	  # }else{
	    # ncp1=(beta*gamma+delta*sigmau)^2*Zadj_sq/
	       # (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  # }
    # }else{
      # print("Wrong input of the sensitivity range.")
      # stop()
    # }

    # if(ncp1<=ncp2){
      # print("Sensitivity range too large")
	  # stop()
    # }else{
      # oldn<-k+2
      # state<-1
      # while(state){
        # temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*oldn)  
        # temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*oldn)  
        # if(temppower < power){
	      # oldn <- oldn*2
	    # }else{
	      # state <- 0
	    # }
      # }  

      # lower <- oldn%/%2
      # upper <- oldn
      # while((upper-lower)>2){
        # new <- (upper+lower)%/%2
        # temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*new)  
        # temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*new)  
        # if(temppower < power){
	      # lower <- new
	    # }else{
	      # upper <- new
	    # }
      # }
    # }
    # return(upper)    
  # }

  # print("Input error.")
  # return(NULL)				  
# }

