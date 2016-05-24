ARsensitivity.power <-function(n,k,lambda,gamma,var.z,sigma1,sigma2,rho,alpha=.05,
                               Delta=NULL, delta=NULL){
if(!is.null(Delta)){
  if(is.numeric(Delta) & length(Delta)==2){
    Delta <- sort(Delta)
    if(!is.null(delta)){
	  if(is.numeric(delta) & length(delta)==1){
	    ncp1 = (lambda*gamma+delta*sigma1)^2*n*var.z/
		       (sigma1^2+2*rho*sigma1*sigma2*lambda+sigma2^2*lambda^2)
		ncp2 = max(Delta^2)*n*var.z
	  }else{
	    print("Wrong input of the true sensitivity parameter")
		stop()
	  }	
	}else{
	  if((Delta[1] < -gamma*lambda/sigma1) & (Delta[2] > -gamma*lambda/sigma1)){
	    ncp1=0
	  }else{
	    ncp1=min((lambda*gamma+Delta[1]*sigma1)^2,(lambda*gamma+Delta[2]*sigma1)^2)*
		     n*var.z/(sigma1^2+2*rho*sigma1*sigma2*lambda+sigma2^2*lambda^2)
	  }
	  ncp2 = max(Delta^2)*n*var.z
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }
}else{
  if(!is.null(delta)){
    print("Sensitivity range missing")
  }else{
    ncp1 = lambda^2*gamma^2*n*var.z/
	       (sigma1^2+2*rho*sigma1*sigma2*lambda+sigma2^2*lambda^2)
	ncp2 = 0
  }
}

temp = qf(1-alpha, df1=1, df2=n-k-1, ncp=ncp2)  
power = 1-pf(temp, df1=1, df2=n-k-1, ncp=ncp1)
return(power)
}
