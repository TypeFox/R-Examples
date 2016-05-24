#'@keywords internal
#'@importFrom stats rgamma runif
sample_scale <- function(c, m, z, U_xi, U_psi,
                         U_Sigma, U_df, ltn, weights, scale){

    n <- length(c)

    #     # Gamma(1/10,2) prior
    #     prior_df <- function(x, log=FALSE){
    #         if(log){
    #             y <- log(x/100)-x/10
    #         }else{
    #             y <- x/100*exp(-x/10)
    #         }
    #         return(y)
    #     }

    #     # Hierarchical prior:
    #     # Gamma(s,2) prior with hyperparameter s ~ Gamma(1,d)
    #     prior_df <- function(x, d, log=FALSE){
    #         if(log){
    #             y <- log(2)+log(d)+log(x)-3*log(x+d)
    #         }else{
    #             y <- 2*d*x/(x+d)^3
    #         }
    #         return(y)
    #     }

    #Jeffrey's prior
    prior_df <- function(x, log=FALSE){
        smalltemp <- trigamma(x/2)-trigamma((x+1)/2)-2*(x+3)/(x*(x+1)^2)
        if(smalltemp <= 0){
            y=ifelse(log,-Inf,0)
        }else{
            if(log){
                y <- 1/2*(log(x)-log(x+3) +log(smalltemp))
            }else{
                y <- sqrt(x/(x+3)*smalltemp)
            }
        }
        return(y)
    }


    c_df <- 2
    fullCl_ind <- which(m!=0)
    fullCl <- length(fullCl_ind)
    acc_rate <- 0

    U_xi_full <- sapply(fullCl_ind, function(j) U_xi[, j])
    U_psi_full <- sapply(fullCl_ind, function(j) U_psi[, j])
    U_Sigma_full <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
    U_df_full <- sapply(fullCl_ind, function(j) U_df[j])

    if(nrow(z)==1){
      U_Sigma_full <- lapply(U_Sigma_full, FUN=matrix, nrow=1, ncol=1)
      U_xi_full <- matrix(U_xi_full, nrow=1)
      U_psi_full <- matrix(U_psi_full, nrow=1)
    }

    df_new <- numeric(n)
    for(j in 1:fullCl){
        df <- U_df_full[j]
        df_new[j] <- 1+exp(stats::runif(1, min=log(df-1)-c_df, max = log(df-1)+c_df))
    }

    loglikold <- mvstlikC(x=z, c=c, clustval=fullCl_ind,
                          xi=U_xi_full, psi=U_psi_full, sigma=U_Sigma_full, df=U_df_full,
                          loglik=TRUE)
    logliknew <- mvstlikC(x=z, c=c, clustval=fullCl_ind,
                          xi=U_xi_full, psi=U_psi_full, sigma=U_Sigma_full, df=df_new,
                          loglik=TRUE)
    u <- stats::runif(fullCl)

    if(fullCl>1){
        for(j in 1:fullCl){


            prob_new <- exp(sum(loglikold$clust[-j],logliknew$clust[j]) + prior_df(df_new[j], log=TRUE) + log(df_new[j]-1)
                            -(loglikold$total + prior_df(U_df_full[j], log=TRUE) + log(U_df_full[j]-1)))
            if(is.na(prob_new)){
                warning("MH probability is NA\n Covariance matrix could be too small, try increasing the hyperprior parameter nu")
                #prob_new <- 0
            }
            if(prob_new==-Inf){
                prob_new <- 0
            }else if(prob_new==Inf){
                prob_new <- 1
            }else{
                prob_new <- min(1, prob_new)
            }
            if (u[j]<prob_new){
                acc_rate <- acc_rate + 1
                U_df_full[j] <- df_new[j]
            }

            obs_j <- which(c==fullCl_ind[j])
            eps <- z[,obs_j, drop=FALSE] - U_xi_full[,j] - sapply(X=ltn[obs_j], FUN=function(x){x*U_psi_full[,j]})
            tra <- traceEpsC(eps, U_Sigma_full[[j]])[,1] # fast C++ code
            #  tra <- apply(X=eps, MARGIN=2, FUN=function(v){sum(diag(tcrossprod(v)%*%solve(U_Sigma_full[[j]])))}) # slow vectorized R code
            scale[obs_j] <- stats::rgamma(length(obs_j),
                                   shape=(U_df[j] + nrow(z) + 1)/2,
                                   rate=(U_df[j] + ltn[obs_j]^2 + tra)/2)
        }
    }else{
        j <- 1
        prob_new <- exp(logliknew$total + prior_df(df_new[j], log=TRUE) + log(df_new[j]-1)
                        -(loglikold$total + prior_df(U_df_full[j], log=TRUE) +log(U_df_full[j]-1)))

        if(prob_new==-Inf){
            warning("MH probability is NA\n Covariance matrix could be too small, try increasing the hyperprior parameter nu")
            #prob_new <- 0
        }else if(prob_new==Inf){
            prob_new <- 1
        }else{
            prob_new <- min(1, prob_new)
        }
        if (u[j]<prob_new){
            acc_rate <- acc_rate + 1
            U_df_full[j] <- df_new[j]
        }

        eps <- z - U_xi_full[,j] - sapply(X=ltn, FUN=function(x){x*U_psi_full[,j]})
        tra <- apply(X=eps, MARGIN=2, FUN=function(v){sum(diag(tcrossprod(v)%*%solve(U_Sigma_full[[j]])))})

        scale <- stats::rgamma(n, shape=(U_df[j] + nrow(z) + 1)/2,
                        rate=(U_df[j] + ltn^2 + tra)/2)
    }


    return(list("df"=U_df_full, "scale"=scale, "acc_rate"=acc_rate/fullCl))
}