##' KFadvanceAR2 function
##'
##' A function to compute one step of the Kalman filter with second order AR state evolution. 
##' Embed in a loop to run the filter on a set of data.
##'
##' The model is: (note that Y and theta are COLUMN VECTORS)
##'
##' theta_t = A*theta_{t-1} + A1*theta_{t-2} + B + C*W (state equation)
##'
##' Y_t = D*theta_t + E + F*V         (observation equation)
##'
##' W and V are the covariance matrices of the state and observation noise. Priors are normal, 
##'
##' N(mu_{t-1},Sigma_{t-1}) and N(mu_{t-2},Sigma_{t-2})
##'
##' Result is the posterior, N(mu_t,Sigma_t), together with the likelihood contribution Prob(Y_t|Y_{t-1})
##'
##' @param obs Y_t
##' @param oldmean mu_{t-1} 
##' @param oldermean mu_{t2} 
##' @param oldvar Sigma_{t-1}
##' @param oldervar Sigma_{t-2}
##' @param A A matrix A
##' @param A1 A matrix A1
##' @param B column vector B
##' @param C matrix C
##' @param D matrix D
##' @param E column vector E
##' @param F matrix F
##' @param W state noise covariance
##' @param V observation noise covariance
##' @param marglik logical, whether to return the marginal likelihood contribution from this observation
##' @param log whether or not to return the log of the likelihood contribution.
##' @param na.rm na.rm logical, whether or not to handle NAs. Defult is FALSE. Set to TRUE if there are any missing values in the observed data.
##' @return list containing the new mean and variance, and if specified, the likelihood
##' @export
KFadvanceAR2 <- function(obs,oldmean,oldermean,oldvar,oldervar,A,A1,B,C,D,E,F,W,V,marglik=FALSE,log=TRUE,na.rm=FALSE){

    if(na.rm){
        if(any(is.na(obs))){
            if(all(is.na(obs))){
                if(log){
                    return(list(mean=A%*%oldmean + A1%*%oldermean  + B,var=A%*%oldvar%*%t(A) + A1%*%oldervar%*%t(A1) + C%*%W%*%t(C),mlik=0))
                }
                else{
                    return(list(mean=A%*%oldmean + A1%*%oldermean  + B,var=A%*%oldvar%*%t(A) + A1%*%oldervar%*%t(A1) + C%*%W%*%t(C),mlik=1))
                }
            }
            else{
                M <- diag(length(obs))
                M <- M[-which(is.na(obs)),]
                obs <- obs[which(!is.na(obs))]
                D <- M%*%D
                E <- M%*%E
                F <- M%*%F
            }
        }
    }
	T <- A%*%oldmean + A1%*%oldermean + B
	S <- A%*%oldvar%*%t(A) + A1%*%oldervar%*%t(A1) + C%*%W%*%t(C)
	K <- D%*%S%*%t(D) + F%*%V%*%t(F)
	if (marglik==TRUE){
		margmean <- D %*% T + E
		if (all(dim(K)==1)){
			newmean <- T + as.numeric(1/K)*S%*%t(D)%*%(obs-margmean)
			newvar <- S - as.numeric(1/K)*S%*%t(D)%*%D%*%S
			marginal <- dnorm(obs,as.numeric(margmean),sqrt(as.numeric(K)),log=log)
		}
		else{
			Kinv <- solve(K)
			newmean <- T + S%*%t(D)%*%Kinv%*%(obs-margmean)
			newvar <- S - S%*%t(D)%*%Kinv%*%D%*%S
			marginal <- dmvnorm(as.vector(obs),as.vector(margmean),K,log=log)
		}
		return(list(mean=newmean,var=newvar,mlik=marginal))
	}
	else{
		if (all(dim(K)==1)){
			newmean <- T + as.numeric(1/K)*S%*%t(D)%*%(obs-D%*%T-E)
			newvar <- S - as.numeric(1/K)*S%*%t(D)%*%D%*%S
		}
		else{
			Kinv <- solve(K)
			newmean <- T + S%*%t(D)%*%Kinv%*%(obs-D%*%T-E)
			newvar <- S - S%*%t(D)%*%Kinv%*%D%*%S
		}
		return(list(mean=newmean,var=newvar))
	}
}
