getQuantilesParams <-
function(fit.i, qFunc=qLNO, quantiles=seq(0.006,0.996,length.out=1000), linksq=c(identity,exp,NULL,NULL),freeParams, fixedParams){
  if(sum(freeParams==c(TRUE, FALSE, FALSE, FALSE))==4){
    params<-linksq[[1]](fit.i$mu.coefficient)
    if(length(fixedParams)==3){
    	samps<-qFunc(quantiles,params[1],fixedParams[1],fixedParams[2],fixedParams[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}
  }else{
    if(sum(freeParams==c(FALSE, TRUE, FALSE, FALSE))==4){
      params<-linksq[[1]](fit.i$sigma.coefficient)
       if(length(fixedParams)==3){
    	samps<-qFunc(quantiles,fixedParams[1],params[1],fixedParams[2],fixedParams[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}
    }else{
      if(sum(freeParams==c(FALSE, FALSE, TRUE, FALSE))==4){
        params<-linksq[[1]](fit.i$nu.coefficient)
        if(length(fixedParams)==3){
    	samps<-qFunc(quantiles,fixedParams[1],fixedParams[2],params[1],fixedParams[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}
    }else{
      if(sum(freeParams==c(FALSE, FALSE, FALSE, TRUE))==4){
        params<-linksq[[1]](fit.i$tau.coefficient)
        if(length(fixedParams)==3){
    	samps<-qFunc(quantiles,fixedParams[1],fixedParams[2],fixedParams[3],params[1])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}
    }else{
      if(sum(freeParams==c(TRUE, TRUE, FALSE, FALSE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$sigma.coefficients))
        if(length(fixedParams)==2){
    	samps<-qFunc(quantiles,params[1],params[2],fixedParams[1],fixedParams[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}
    }else{
      if(sum(freeParams==c(TRUE, FALSE, TRUE, FALSE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$nu.coefficients))
        if(length(fixedParams)==2){
    	samps<-qFunc(quantiles,params[1],fixedParams[1],params[2],fixedParams[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(TRUE, FALSE, FALSE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$tau.coefficients))
        if(length(fixedParams)==2){
        	    	samps<-qFunc(quantiles,params[1],fixedParams[1],fixedParams[2],params[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}  
    }else{
      if(sum(freeParams==c(FALSE, TRUE, TRUE, FALSE))==4){
        params<-c(linksq[[1]](fit.i$sigma.coefficients), linksq[[2]](fit.i$nu.coefficients))
        if(length(fixedParams)==2){
    	samps<-qFunc(quantiles,fixedParams[1],params[1],params[2],fixedParams[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}  
    }else{
      if(sum(freeParams==c(FALSE, TRUE, FALSE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$sigma.coefficients), linksq[[2]](fit.i$tau.coefficients))
        if(length(fixedParams)==2){
    	samps<-qFunc(quantiles,fixedParams[1],params[1],fixedParams[2],params[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}}  
    }else{
      if(sum(freeParams==c(FALSE, FALSE, TRUE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$nu.coefficients), linksq[[2]](fit.i$tau.coefficients))
        if(length(fixedParams)==2){
    	samps<-qFunc(quantiles,fixedParams[1],fixedParams[2],params[1],params[2])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(TRUE, TRUE, TRUE, FALSE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$sigma.coefficients), linksq[[3]](fit.i$nu.coefficients))
        if(length(fixedParams)==1){
    	samps<-qFunc(quantiles,params[1],params[2],params[3],fixedParams[1])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2],params[3])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(TRUE, TRUE, FALSE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$sigma.coefficients), linksq[[3]](fit.i$tau.coefficients))
        if(length(fixedParams)==1){
    	samps<-qFunc(quantiles,params[1],params[2],fixedParams[1],params[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2],params[3])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(TRUE, FALSE, TRUE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$nu.coefficients), linksq[[3]](fit.i$tau.coefficients))
        if(length(fixedParams)==1){
    	samps<-qFunc(quantiles,params[1],fixedParams[1],params[2],params[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2],params[3])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(FALSE, TRUE, TRUE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$sigma.coefficients), linksq[[2]](fit.i$nu.coefficients), linksq[[3]](fit.i$tau.coefficients))
        if(length(fixedParams)==1){
    	samps<-qFunc(quantiles,fixedParams[1],params[1],params[2],params[3])
    	}else{
    		if(length(fixedParams)==0){
    			samps<-qFunc(quantiles,params[1],params[2],params[3])
    		}else{
    			samps<-rep(NA, length(quantiles))
    			cat('problem with fixed and free parameters, quantiles are wrong','\n','\n')
    			}} 
    }else{
      if(sum(freeParams==c(TRUE, TRUE, TRUE, TRUE))==4){
        params<-c(linksq[[1]](fit.i$mu.coefficients), linksq[[2]](fit.i$sigma.coefficients), linksq[[3]](fit.i$nu.coefficients),linksq[[4]](fit.i$tau.coefficients))
        samps<-qFunc(quantiles,params[1],params[2],params[3],params[4])
    }else{
      warning('There is a problem, you have',length(fit.i$parameters),'parameters and that is either too few or too many!  If it is more than 4, contact the function author.','\n','\n')
      samps<-NA
      params<-NA
    }#end final else
    }}}}}}}}}}}}}}#end previous 14 elses
  out<-list('samps'=samps,'params'=params)
  return(out)
}
