#'EM MAP for mixture of sNiW
#'
#'Maximum A Posteriori (MAP) estimation of mixture of
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@details \code{MAP_sNiW_mmEM} provides an estimation for the MAP of mixtures of
#'Normal inverse Wishart distributed observations. \code{MAP_sNiW_mmEM_vague} provides
#'an estimates incorporating a vague component in the mixture.
#'\code{MAP_sNiW_mmEM_weighted} provides a weigthed version of the algorithm.
#'
#'@param xi_list a list of length \code{n}, each element is a vector of size \code{d}
#'containing the argument \code{xi} of the corresponding allocated cluster.
#'
#'@param psi_list a list of length \code{n}, each element is a vector of size \code{d}
#'containing the argument \code{psi} of the corresponding allocated cluster.
#'
#'@param S_list a list of length \code{n}, each element is a matrix of size \code{d x d}
#'containing the argument \code{S} of the corresponding allocated cluster.
#'
#'@param hyperG0 prior mixing distribution used if \code{init} is \code{NULL}.
#'
#'@param init a list for initializing the algorithm with the following elements: \code{b_xi},
#'\code{b_psi}, \code{lambda}, \code{B}, \code{nu}. Default is \code{NULL} in which case
#'the initialization of the algorithm is random.
#'
#'@param K integer giving the number of mixture components.
#'
#'@param maxit integer giving the maximum number of iteration for the EM algorithm.
#'Default is \code{100}.
#'
#'@param tol real number giving the tolerance for the stopping of the EM algorithm.
#'Default is \code{0.1}.
#'
#'@param doPlot a logical flag indicating wether the algorithm progression should be plotted.
#'Default is \code{TRUE}.
#'
#'@param verbose logical flag indicating whether plot should be drawn. Default is \code{TRUE}.
#'
#'@rdname MAP_sNiW_mmEM
#'
#'@author Boris Hejblum, Chariff Alkhassim
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
#'
#'@export
MAP_sNiW_mmEM<- function(xi_list, psi_list, S_list, hyperG0, init=NULL, K, maxit=100, tol=1E-1,
                          doPlot=TRUE, verbose=TRUE){



  N <- length(xi_list)
  d <- length(xi_list[[1]])

  if(length(psi_list) != N | length(S_list) != N){
    stop("Number of observations/MCMC iterations is not matching")
  }

  U_xi <- list() #matrix(0, nrow=d,ncol=K)
  U_psi <- list() #matrix(0, nrow=d,ncol=K)
  U_Sigma <- list() # array(dim=c(d,d,K))
  U_B <- list() #array(dim=c(2,2,K))
  U_df <- list() #numeric(K)


  #priors
  #   alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  xi_p <- apply(sapply(xi_list, "["), MARGIN=1, FUN=mean)
  psi_p <- apply(sapply(psi_list, "["), MARGIN=1, FUN=mean)
  kappa0 <- 0.01
  nu<- d+1
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation


  if(is.null(init)){
    for(k in 1:K){
      #sampling the cluster parameters
      NNiW <- rNNiW(hyperG0, diagVar=FALSE)
      U_xi[[k]] <- NNiW[["xi"]]
      U_psi[[k]] <- NNiW[["psi"]]
      U_Sigma[[k]] <- NNiW[["S"]]
      U_B[[k]] <- diag(100, 2)
      U_df[[k]] <- d+1
    }
  }
  else{
    for(k in 1:K){
      #cluster parameters
      U_xi[[k]] <- init[[k]]$b_xi
      U_psi[[k]] <- init[[k]]$b_psi
      U_Sigma[[k]] <- init[[k]]$lambda
      U_B[[k]] <- init[[k]]$B
      U_df[[k]] <- init[[k]]$nu
    }
  }


  K<-K
  #priors
  alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  weights <- rep(1/K, K)

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


  for(i in 1:maxit){

    #E step

    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))


    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


    #M step
    N_k <- rowSums(r)
    weights  <- (N_k + alpha[k] - 1)/(N + sum(alpha) - K)

    for(k in 1:K){


      Sinv_list <- mapply(S = S_list,
                          FUN=function(S){solve(S)},
                          SIMPLIFY=FALSE)
      Sinv_sum <- Reduce('+', Sinv_list)
      rSinv_list <- mapply(S = Sinv_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*S},
                           SIMPLIFY=FALSE)
      rSinv_sum <- Reduce('+', rSinv_list)
      r_xi<-apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_psi<-apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_xi_psi<-cbind(r_xi,r_psi)
      U_xi_U_psi<-((kappa0/N*t(matrix(cbind(xi_p,psi_p),d,d))%*%Sinv_sum)+Reduce('+',lapply(1:N,
                                                                                            function(i,m) {t(matrix(r_xi_psi[i,],d,d))%*%m[[i]]},
                                                                                            m=rSinv_list)))%*% solve(kappa0/N*Sinv_sum+rSinv_sum)

      U_xi[[k]] <- U_xi_U_psi[1,]
      U_psi[[k]]<- U_xi_U_psi[2,]

      xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
      psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
      xim0 <- U_xi[[k]] - xi_p
      psim0 <- U_xi[[k]] - psi_p

      Sinv_list <- lapply(S_list,solve)
      Sinv_sum <- Reduce('+', Sinv_list)
      rSinv_list <- mapply(Sinv = Sinv_list,
                           rik = as.list(r[k, ]),
                           FUN=function(Sinv, rik){rik*Sinv},
                           SIMPLIFY=FALSE)
      rSinv_sum <- Reduce('+', rSinv_list)
      U_B[[k]] <- 1/(N_k[k]*d + d + 1)*(solve(C) + matrix(rowSums(mapply(x = xim,  #1/
                                                                         p = psim,
                                                                         rSinv = rSinv_list,
                                                                         FUN=function(x,p,rSinv){
                                                                           v <- rbind(x, p)
                                                                           v%*%rSinv%*%t(v)
                                                                         }, SIMPLIFY=TRUE)),
                                                          nrow=2, byrow=FALSE)
                                        +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
      )
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)

      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                     - N_k[k]*d*log(N_k[k]*nu0/2)
                                                     + const_nu0_uniroot
      )}, lower = d, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}



      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }
    #        cat("df",unlist(U_df), "\n")

    temp_logliks <- log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                     U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                     U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE)))
    loglik[i+1] <- sum(temp_logliks)


    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      #       browser()
      temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
      loglik[i+1] <- sum(temp_logliks)
    }

    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }

    if(doPlot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
                     ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }

    if(abs(loglik[i+1]-loglik[i])<tol){break}
  }
  U_B <- lapply(U_B,function(i)(diag(diag(i))))
  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}


#'@rdname MAP_sNiW_mmEM
#'
#'@param obsweight_list a list of length \code{n} where each element is a vector of weights for
#'each sampled cluster at each MCMC iterations.
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
#'
#'@export
MAP_sNiW_mmEM_weighted<- function(xi_list, psi_list, S_list, obsweight_list,
                                   hyperG0, K, maxit=100, tol=1E-1, doPlot=TRUE, verbose=TRUE){


  pseudoN <- length(xi_list)
  N <- sum(rep(1,pseudoN)*unlist(obsweight_list))
  d <- length(hyperG0[[1]])

  if(length(psi_list) != pseudoN | length(S_list) != pseudoN | length(obsweight_list) != pseudoN){
    stop("Number of observations/MCMC iterations is not matching")
  }

  U_xi <- list() #matrix(0, nrow=d,ncol=K)
  U_psi <- list() #matrix(0, nrow=d,ncol=K)
  U_Sigma <- list() # array(dim=c(d,d,K))
  U_B <- list() #array(dim=c(2,2,K))
  U_df <- list() #numeric(K)

  #priors
  alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  xi_p <- apply(sapply(xi_list, "["), MARGIN=1, FUN=mean)
  psi_p <- apply(sapply(psi_list, "["), MARGIN=1, FUN=mean)
  kappa0 <- 0.01
  nu<- d+1
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters
    NNiW <- rNNiW(hyperG0, diagVar=FALSE)
    U_xi[[k]] <- NNiW[["xi"]]
    U_psi[[k]] <- NNiW[["psi"]]
    U_Sigma[[k]] <- NNiW[["S"]]
    U_B[[k]] <- diag(100, 2)
    U_df[[k]] <- d+1
  }

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))

  for(i in 1:maxit){

    #E step
    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


    #M step
    N_k <- rowSums(r)
    weights  <- (N_k + alpha[k] - 1)/(N + sum(alpha) - K)

    for(k in 1:K){
      xi_m_k_xNk <- colSums(apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
      U_xi[[k]] <- (xi_m_k_xNk + kappa0/N*xi_p)/(N_k[k]  + kappa0/N)
      psi_m_k_xNk <- colSums(apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))
      U_psi[[k]] <- (psi_m_k_xNk + kappa0/N*psi_p)/(N_k[k] + kappa0/N)

      xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
      psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
      xim0 <- U_xi[[k]] - xi_p
      psim0 <- U_xi[[k]] - psi_p

      Sinv_list <- lapply(S_list,solve)
      Sinv_sum <- Reduce('+', Sinv_list)
      rSinv_list <- mapply(Sinv = Sinv_list,
                           rik = as.list(r[k, ]),
                           FUN=function(Sinv, rik){rik*Sinv},
                           SIMPLIFY=FALSE)
      rSinv_sum <- Reduce('+', rSinv_list)
      U_B[[k]] <- 1/(N_k[k]*d + d + 1)*solve(solve(C) + matrix(rowSums(mapply(x = xim,
                                                                              p = psim,
                                                                              rSinv = rSinv_list,
                                                                              FUN=function(x,p,rSinv){
                                                                                v <- rbind(x, p)
                                                                                v%*%rSinv%*%t(v)
                                                                              }, SIMPLIFY=TRUE)),
                                                               nrow=2, byrow=FALSE)
                                             +kappa0/N*rbind(xim0, psim0)%*%Sinv_sum%*%t(rbind(xim0, psim0))
      )
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)
      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                     - N_k[k]*d*log(N_k[k]*nu0/2)
                                                     + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}


      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }

    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      #browser()
    }
    if(abs(loglik[i+1]-loglik[i])<tol){break}

    if(doPlot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
                     ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }

  if(doPlot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
                   ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
  }

  return(list("r"=r,
              "loglik" = loglik[2:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}









#'@rdname MAP_sNiW_mmEM
#'
#'@importFrom stats var uniroot
#'
#'@importFrom graphics plot
#'
#'@export
#'
#'@examples
#'set.seed(1234)
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(0.3, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 100
#'hyperG0$D_psi <- 100
#'hyperG0$nu <- 20
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(1, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.1
#'hyperG0$D_xi <- 1
#'hyperG0$D_psi <- 1
#'hyperG0$nu <- 2
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'w_list <- list()
#'
#'for(k in 1:200){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'  w_list [[k]] <- 0.75
#'}
#'
#'
#'hyperG02 <- list()
#'hyperG02$b_xi <- c(-1, 2)
#'hyperG02$b_psi <- c(-0.1, 0.5)
#'hyperG02$kappa <- 0.1
#'hyperG02$D_xi <- 1
#'hyperG02$D_psi <- 1
#'hyperG02$nu <- 4
#'hyperG02$lambda <- 0.5*diag(2)
#'
#'for(k in 201:400){
#'  NNiW <- rNNiW(hyperG02, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'  w_list [[k]] <- 0.25
#'
#'}
#'
#'map <- MAP_sNiW_mmEM(xi_list, psi_list, S_list, hyperG0, K=2, tol=0.1)
#'
MAP_sNiW_mmEM_vague <- function(xi_list, psi_list, S_list,
                                 hyperG0, K=10, maxit=100, tol=1E-1, doPlot=TRUE, verbose=TRUE){

  N <- length(xi_list)
  d <- length(hyperG0[[1]])

  if(length(psi_list) != N | length(S_list) != N){
    stop("Number of observations/MCMC iterations is not matching")
  }

  U_xi <- list() #matrix(0, nrow=d,ncol=K)
  U_psi <- list() #matrix(0, nrow=d,ncol=K)
  U_Sigma <- list() # array(dim=c(d,d,K))
  U_B <- list() #array(dim=c(2,2,K))
  U_df <- list() #numeric(K)


  #priors
  alpha <- rep(1, K) #parameters of a Dirichlet prior on the cluster weights
  nu<- d+1
  lambda<- diag(apply(sapply(xi_list, "["),MARGIN=1, FUN=stats::var))
  C <- diag(2)*1000
  L <- (diag(apply(sapply(xi_list, "["), MARGIN=1, FUN=stats::var))
        + diag(apply(sapply(psi_list, "["), MARGIN=1, FUN=stats::var))
  )/2


  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters
    NNiW <- rNNiW(hyperG0, diagVar=FALSE)
    U_xi[[k]] <- NNiW[["xi"]]
    U_psi[[k]] <- NNiW[["psi"]]
    U_Sigma[[k]] <- NNiW[["S"]]
    U_B[[k]] <- diag(100, 2)
    U_df[[k]] <- d+1
  }


  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))
  #Q <- numeric(maxit+1)
  #Q[1] <- -Inf

  for(i in 1:maxit){
    #browser()

    #E step
    r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))



    #M step
    N_k <- rowSums(r)
    weights  <- N_k/N #(N_k + alpha[k] - 1)/(N + sum(alpha) - K) #equivalent for alpha[k]=1

    for(k in 1:K){
      U_xi[[k]] <- colSums(apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))/N_k[k]
      U_psi[[k]] <- colSums(apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x}))/N_k[k]

      xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
      psim <- lapply(psi_list, function(x){x - U_psi[[k]]})
      rSinv_list <- mapply(S = S_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*solve(S)},
                           SIMPLIFY=FALSE)
      rSinv_sum <- Reduce('+', rSinv_list)
      U_B[[k]] <- 1/(N_k[k]*d + 1)*(solve(C) + matrix(rowSums(mapply(x = xim,
                                                                     p = psim,
                                                                     rSinv = rSinv_list,
                                                                     FUN=function(x,p,rSinv){
                                                                       v <- rbind(x, p)
                                                                       v%*%rSinv%*%t(v)
                                                                     }, SIMPLIFY=TRUE)),
                                                      nrow=2, byrow=FALSE))
      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)
      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                     - N_k[k]*d*log(N_k[k]*nu0/2)
                                                     + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}


      U_Sigma[[k]] <- (N_k[k]*U_df[[k]] + 1)*solve(solve(L) + rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))



    if(is.na(loglik[i+1]) | is.nan(loglik[i+1]) | is.infinite(loglik[i+1])){
      temp_logliks[which(is.infinite(temp_logliks))] <- min(temp_logliks[-which(is.infinite(temp_logliks))])
      loglik[i+1] <- sum(temp_logliks)
    }

    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }
    if(abs(loglik[i+1]-loglik[i])<tol){break}

    if(doPlot){
      graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
                     ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
    }
  }

  if(doPlot){
    graphics::plot(y=loglik[2:(i+1)], x=c(1:i),
                   ylab="Log-likelihood", xlab="Iteration", type="b", col="blue", pch=16)
  }

  return(list("r"=r,
              "loglik" = loglik[1:(i+1)],
              "U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma,
              "weights"=weights))

}
