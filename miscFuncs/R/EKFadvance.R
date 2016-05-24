##' EKFadvance function
##'
##' A function to perform one iteration of ther EKF. Currently UNDER DEVELOPMENT.
##'
##' @param obs  observations
##' @param oldmean old mean
##' @param oldvar old variance
##' @param phi Function computing a Taylor Series approximation of the system equation. Can include higher (ie 2nd order and above) terms.
##' @param phi.arglist arguments for function phi
##' @param psi Function computing a Taylor Series approximation of the observation equation. Can include higher (ie 2nd order and above) terms.
##' @param psi.arglist arguments for function psi
##' @param W system noise matrix
##' @param V observation noise matrix
##' @param loglik whether or not to compute the pseudo-likelihood
##' @param na.rm logical, whether or not to handle NAs. Defult is FALSE. Set to TRUE if there are any missing values in the observed data.
##' @return list containing the new mean and variance, and if specified, the likelihood
##' @export
EKFadvance <- function(obs,oldmean,oldvar,phi,phi.arglist,psi,psi.arglist,W,V,loglik=FALSE,na.rm=FALSE){ # phi==phi(Theta_t-1,W_t) but specified as phi(Theta_t-1) since only ever evaluate phi(Theta_t-1,0) ... similar for psi 

    s <- phi(oldmean,phi.arglist) # s = system equation bits
    if(is.null(s$higherorder)){
        s$higherorder <- list(mean=0,var=0)
    }
    S <- s$dtheta%*%oldvar%*%t(s$dtheta) + s$dW%*%W%*%t(s$dW) + s$higherorder$var # Sigma_t|t-1 WITH higher order terms    
    
    o <- psi(s$funval + s$higherorder$mean,psi.arglist,S=S,s=s) # o = observation equation bits WITH higher order terms
    if(is.null(o$higherorder)){
        o$higherorder <- list(mean=0,var=0)
    }  
        
    if(na.rm){
        if(any(is.na(obs))){
            if(all(is.na(obs))){
                return(list(mean=s$funval + s$higherorder$mean,var=S,loglik=0))
            }
            else{
                M <- diag(length(obs))
                M <- M[-which(is.na(obs)),]
                obs <- obs[which(!is.na(obs))]
                o$funval <- M%*%o$funval
                o$dtheta <- M%*%o$dtheta
                o$dV <- M%*%o$dV
                if(!identical(o$higherorder$mean,0)){ # leads to a small computational saving if higher order terms are not used in the model
                    o$higherorder$mean <- M%*%o$higherorder$mean
                }
                if(!identical(o$higherorder$mean,0)){                
                    o$higherorder$var <- M%*%o$higherorder$var%*%t(M)
                }
            }
        }
    }
    
    K <- as.matrix(o$dtheta%*%S%*%t(o$dtheta) + o$dV%*%V%*%t(o$dV) + o$higherorder$var)
            
	Kinv <- solve(K)
	newmean <- s$funval + S%*%t(o$dtheta)%*%Kinv%*%(obs-o$funval-o$higherorder$mean)
    newvar <- S - S%*%t(o$dtheta)%*%Kinv%*%o$dtheta%*%S
    if(loglik){  
	    ll <- dmvnorm(as.vector(obs),as.vector(o$funval),K,log=TRUE)
	}
		
	retlist <- list(mean=newmean,var=newvar)
	if(loglik){
	    retlist$loglik <- ll
	}	
	return(retlist)
}
