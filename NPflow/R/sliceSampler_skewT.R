#'@keywords internal
#'@author Boris Hejblum
#'@importFrom stats rbeta rgamma runif
sliceSampler_skewT <- function(c, m, alpha, z, hyperG0, U_xi, U_psi,
                               U_Sigma, U_df, scale, diagVar){


    maxCl <- length(m) #maximum number of clusters
    ind <- which(m!=0) # indexes of non empty clusters

    # Sample the weights, i.e. the frequency of each existing cluster from a Dirichlet:
    # temp_1 ~ Gamma(m_1,1), ... , temp_K ~ Gamma(m_K,1)    # and sample the rest of the weigth for potential new clusters:
    # temp_{K+1} ~ Gamma(alpha, 1)
    # then renormalise temp
    w <- numeric(maxCl)
    temp <- stats::rgamma(n=(length(ind)+1), shape=c(m[ind], alpha), scale = 1)
    temp_norm <- temp/sum(temp)
    w[ind] <- temp_norm[-length(temp_norm)]
    R <- temp_norm[length(temp_norm)]
    #R is the rest, i.e. the weight for potential new clusters

    # Sample the latent u
    u  <-stats::runif(maxCl)*w[c]
    u_star <- min(u)

    # Sample the remaining weights that are needed with stick-breaking
    # i.e. the new clusters
    ind_new <- which(m==0) # potential new clusters
    if(length(ind_new)>0){
        t <- 0 # the number of new non empty clusters
        while(R>u_star && (t<length(ind_new))){
            # sum(w)<1-min(u) <=> R>min(u) car R=1-sum(w)
            t <- t+1
            beta_temp <-stats::rbeta(n=1, shape1=1, shape2=alpha)
            # weight of the new cluster
            w[ind_new[t]] <- R*beta_temp
            R <- R * (1-beta_temp) # remaining weight
        }
        ind_new <- ind_new[1:t]

        # Sample the centers and spread of each new cluster from prior
        for (i in 1:t){
            NNiW <- rNNiW(hyperG0, diagVar)
            U_xi[, ind_new[i]] <- NNiW[["xi"]]
            U_psi[, ind_new[i]] <- NNiW[["psi"]]
            U_Sigma[, , ind_new[i]] <- NNiW[["S"]]
            U_df[ind_new[i]] <- 10
        }
    }

    fullCl_ind <- which(w != 0)
    # likelihood of belonging to each cluster computation
    # sampling clusters
    U_xi_full <- sapply(fullCl_ind, function(j) U_xi[, j])
    U_psi_full <- sapply(fullCl_ind, function(j) U_psi[, j])
    U_Sigma_full <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
    U_df_full <- sapply(fullCl_ind, function(j) U_df[j])

    if(nrow(z)==1){
      U_Sigma_full <- lapply(U_Sigma_full, FUN=matrix, nrow=1, ncol=1)
      U_xi_full <- matrix(U_xi_full, nrow=1)
      U_psi_full <- matrix(U_psi_full, nrow=1)
    }

    if(length(fullCl_ind)>1){
        l <- mmvstpdfC(x=z, xi=U_xi_full, psi=U_psi_full, sigma=U_Sigma_full, df=U_df_full, Log = FALSE)
        u_mat <- t(sapply(w[fullCl_ind], function(x){as.numeric(u < x)}))
        prob_mat <- u_mat * l

        #fast C++ code
        c <- fullCl_ind[sampleClassC(prob_mat)]
        #         #slow C++ code
        #         c <- fullCl_ind[sampleClassC_bis(prob_mat)]
        #         #vectorized R code
        #         c <- fullCl_ind[apply(X= prob_mat, MARGIN=2, FUN=function(v){match(1,rmultinom(n=1, size=1, prob=v))})]
        #         #alternative implementation:
        #         prob_colsum <- colSums(prob_mat)
        #         prob_norm <- apply(X=prob_mat, MARGIN=1, FUN=function(r){r/prob_colsum})
        #         c <- fullCl_ind[apply(X=prob_norm, MARGIN=1, FUN=function(r){match(TRUE,runif(1) <cumsum(r))})]
    }else{
        c <- rep(fullCl_ind, maxCl)
    }

    #vectorized R code
    #     U_xi_list <- lapply(fullCl_ind, function(j) U_xi[, j])
    #     U_psi_list <- lapply(fullCl_ind, function(j) U_psi[, j])
    #     U_Sigma_list <- lapply(fullCl_ind, function(j) U_Sigma[, ,j])
    #     U_df_list <- lapply(fullCl_ind, function(j) U_df[j])
    #     l <- mvstpdf(x=z, xi=U_xi_list, sigma=U_Sigma_list, psi=U_psi_list, df=U_df_list)

    # non vectorized code for cluster allocation:
    #     nb_fullCl_ind <- length(fullCl_ind)
    #     l <- numeric(nb_fullCl_ind) # likelihood of belonging to each cluster
    #     for(i in 1:maxCl){
    #         for (j in fullCl_ind){
    #             l[j] <- mvsnpdf(x = matrix(z[,i], ncol= 1, nrow=length(z[,i])) ,
    #                             xi = U_xi[, j],
    #                             sigma = U_Sigma[, , j],
    #                             psi = U_psi[,j]
    #
    #             )
    #         }
    #     }

    m_new <- numeric(maxCl) # number of observations in each cluster
    m_new[unique(c)] <- table(c)[as.character(unique(c))]

    ltn <- numeric(maxCl) # latent truncated normal variables
    for (k in which(m_new!=0)){
        obs_k <- which(c==k)
        siginv <- solve(U_Sigma[, , k])
        psi <- U_psi[,k, drop=FALSE]
        A_k <-  1/(1 + (crossprod(psi, siginv)%*%psi))
        a_ik <- (tcrossprod(A_k, psi)%*%siginv%*%(z[,obs_k, drop=FALSE]-U_xi[,k]))
        A_k <- A_k/scale[obs_k]
        ltn[obs_k] <- rtruncnorm(length(obs_k), a=0, b=Inf, mean = a_ik, sd = sqrt(A_k))
    }

    return(list("c"=c, "m"=m_new, "weights"=w, "latentTrunc"=ltn,
                "xi"=U_xi, "psi"=U_psi, "Sigma"=U_Sigma, "df"=U_df))
}