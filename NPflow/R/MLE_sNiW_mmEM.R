#'EM MLE for mixture of sNiW
#'
#'Maximum likelihood estimation of mixture of
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@param xi_list a list of length \code{N} whose elements are observed vectors of length \code{d}
#'of the mean parameters xi.
#'
#'@param psi_list a list of length \code{N} whose elements are observed vectors of length \code{d}
#'of the skew parameters psi.
#'
#'@param S_list a list of length \code{N} whose elements are observed variance-covariance matrices
#'of dimension \code{d x d}.
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
#'@param doPlot a logical flag indicating wether the algorithm progression should be plotted. Default is \code{TRUE}.
#'
#'@param verbose logical flag indicating whether plot should be drawn. Default is \code{TRUE}.
#'
#'@rdname MLE_sNiW_mmEM
#'
#'@importFrom stats uniroot
#'
#'@importFrom graphics plot
#'
#'@export MLE_sNiW_mmEM
#'
#'@author Boris Hejblum, Chariff Alkhassim
#'
#'@examples
#'set.seed(1234)
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(0.3, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 100
#'hyperG0$D_psi <- 100
#'hyperG0$nu <- 3
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:200){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'hyperG02 <- list()
#'hyperG02$b_xi <- c(-1, 2)
#'hyperG02$b_psi <- c(-0.1, 0.5)
#'hyperG02$kappa <- 0.001
#'hyperG02$D_xi <- 10
#'hyperG02$D_psi <- 10
#'hyperG02$nu <- 3
#'hyperG02$lambda <- 0.5*diag(2)
#'
#'for(k in 201:400){
#'  NNiW <- rNNiW(hyperG02, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'mle <- MLE_sNiW_mmEM(xi_list, psi_list, S_list, hyperG0, K=2)
#'
MLE_sNiW_mmEM <- function(xi_list, psi_list, S_list, hyperG0, K, init=NULL, maxit=100, tol=1E-1,
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
      U_Sigma[[k]] <-  init[[k]]$lambda
      U_B[[k]] <- init[[k]]$B
      U_df[[k]] <- init[[k]]$nu
    }
  }

  weights <- rep(1/K, K)
  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                    U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                    U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))



  #Q <- numeric(maxit+1)
  #Q[1] <- -Inf


  r <- mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                  U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                  U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["))

  logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
  r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
  rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
  r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))
  N_k<-rowSums(r)




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
    weights  <- N_k/N
    #     cat("weights:", weights, "\n")


    for(k in 1:K){

      rSinv_list <- mapply(S = S_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*solve(S)},
                           SIMPLIFY=FALSE)

      rSinv_sum <- Reduce('+', rSinv_list)


      r_xi<-apply(X=sapply(xi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})
      r_psi<-apply(X=sapply(psi_list, FUN="["), MARGIN=1, FUN=function(x){r[k, ]*x})

      r_xi_psi<-cbind(r_xi,r_psi)


      U_xi_U_psi<-Reduce('+',lapply(1:N,
                                    function(i,m) {t(matrix(r_xi_psi[i,],d,d))%*%m[[i]]},
                                    m=rSinv_list))%*% solve(rSinv_sum)


      U_xi[[k]] <- U_xi_U_psi[1,]
      U_psi[[k]]<- U_xi_U_psi[2,]
      xim <- lapply(xi_list, function(x){x - U_xi[[k]]})
      psim <- lapply(psi_list, function(x){x - U_psi[[k]]})

      #       U_B[[k]] <- 1/(N_k[k]*d)*(matrix(rowSums(mapply(x = xim,
      #                                                       p = psim,
      #                                                       rSinv = rSinv_list,
      #                                                       FUN=function(x,p,rSinv){
      #                                                         v <- rbind(x, p)
      #                                                         v%*%rSinv%*%t(v)
      #                                                       }, SIMPLIFY=TRUE)),
      #                                        nrow=2, byrow=FALSE))

      U_B[[k]] <- 1/(N_k[k]*d)*(diag(diag((matrix(rowSums(mapply(x = xim,
                                                                 p = psim,
                                                                 rSinv = rSinv_list,
                                                                 FUN=function(x,p,rSinv){
                                                                   v <- rbind(x, p)
                                                                   v%*%rSinv%*%t(v)
                                                                 }, SIMPLIFY=TRUE)),
                                                  nrow=2, byrow=FALSE)))))




      const_nu0_uniroot <- (sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                            + N_k[k]*log(det(rSinv_sum))
                            + 2)

      U_df[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                     - N_k[k]*d*log(N_k[k]*nu0/2)
                                                     + const_nu0_uniroot
      )}, lower = d+1, upper=1E12)$root, TRUE)
      if(inherits(U_df[[k]], "try-error")){U_df[[k]] <- d+1}

      U_Sigma[[k]] <- N_k[k]*U_df[[k]]*solve(rSinv_sum)
    }

    loglik[i+1] <- sum(log(colSums(weights%*%mmsNiWpdfC(xi = sapply(xi_list, "["), psi = sapply(psi_list, "["), Sigma = S_list,
                                                        U_xi0 = sapply(U_xi, "["), U_psi0 = sapply(U_psi, "["), U_B0 =U_B,
                                                        U_Sigma0 = U_Sigma, U_df0 = sapply(U_df, "["), Log=FALSE))))


    #     cat("it ", i, ": loglik = ", loglik[i+1],"\n\n", sep="")
    if(abs(loglik[i+1]-loglik[i])<tol){break}
    if(verbose){
      cat("it ", i, ": loglik = ", loglik[i+1],"\n", sep="")
      cat("weights:", weights, "\n\n")
    }

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




#'EM MLE for mixture of NiW
#'
#'Maximum likelihood estimation of mixture of
#'Normal inverse Wishart distributed observations with an EM algorithm
#'
#'@param mu_list a list of length \code{N} whose elements are observed vectors of length \code{d}
#'of the mean parameters.
#'
#'@param S_list a list of length \code{N} whose elements are observed variance-covariance matrices
#'of dimension \code{d x d}.
#'
#'@param hyperG0 prior mixing distribution used for randomly initializing the algorithm.
#'
#'@param K integer giving the number of mixture components.
#'
#'@param maxit integer giving the maximum number of iteration for the EM algorithm.
#'Default is \code{100}.
#'
#'@param tol real number giving the tolerance for the stopping of the EM algorithm.
#'Default is \code{0.1}.
#'
#'@param doPlot a logical flag indicating wether the algorithm progression should be plotted. Default is \code{TRUE}.
#'
#'@rdname MLE_NiW_mmEM
#'
#'@importFrom stats uniroot
#'
#'@importFrom graphics plot
#'
#'@export
#'
#'@examples
#'set.seed(123)
#' U_mu <- list()
#' U_Sigma <- list()
#' U_nu<-list()
#' U_kappa<-list()
#'
#' d <- 2
#' hyperG0 <- list()
#' hyperG0[["mu"]] <- rep(1,d)
#' hyperG0[["kappa"]] <- 0.01
#' hyperG0[["nu"]] <- d+1
#' hyperG0[["lambda"]] <- diag(d)
#'
#' for(k in 1:200){
#'
#'   NiW <- rNiW(hyperG0, diagVar=FALSE)
#'   U_mu[[k]] <-NiW[["mu"]]
#'   U_Sigma[[k]] <-NiW[["S"]]
#' }
#'
#'
#' hyperG02 <- list()
#' hyperG02[["mu"]] <- rep(2,d)
#' hyperG02[["kappa"]] <- 1
#' hyperG02[["nu"]] <- d+10
#' hyperG02[["lambda"]] <- diag(d)/10
#'
#' for(k in 201:400){
#'
#'   NiW <- rNiW(hyperG02, diagVar=FALSE)
#'   U_mu[[k]] <-NiW[["mu"]]
#'   U_Sigma[[k]] <-NiW[["S"]]
#' }
#'
#'
#' mle <- MLE_NiW_mmEM( U_mu, U_Sigma, hyperG0, K=2)
#'
#' hyperG0[["mu"]]
#' hyperG02[["mu"]]
#' mle$U_mu
#'
#' hyperG0[["lambda"]]
#' hyperG02[["lambda"]]
#' mle$U_lambda
#'
#' hyperG0[["nu"]]
#' hyperG02[["nu"]]
#' mle$U_nu
#'
#' hyperG0[["kappa"]]
#' hyperG02[["kappa"]]
#' mle$U_kappa


MLE_NiW_mmEM <- function(mu_list, S_list, hyperG0, K, maxit=100, tol=1e-1, doPlot=TRUE){

  N <- length(mu_list)
  d <- length(hyperG0[[1]])

  if(length(mu_list) != N | length(S_list) != N){
    stop("Number of observations/MCMC iterations is not matching")
  }

  U_mu <- list()
  U_Sigma <- list()
  U_nu<-list()
  U_kappa<-list()

  #initialisation
  weights <- rep(1/K, K)
  for(k in 1:K){
    #sampling the cluster parameters

    NiW <- rNiW(hyperG0, diagVar=FALSE)
    U_mu[[k]] <-NiW[["mu"]]
    U_Sigma[[k]] <-NiW[["S"]]
    U_nu[[k]]<-d+1
    U_kappa[[k]]<-0.01
  }

  loglik <- numeric(maxit+1)
  loglik[1] <- sum(log(colSums(weights%*%mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
                                                  U_mu0=sapply(U_mu, "["),
                                                  U_kappa0=sapply(U_kappa, "["),
                                                  U_nu0=sapply(U_nu, "["),
                                                  U_lambda0=U_Sigma,Log=FALSE))))
  loglik_w <- loglik

  for(i in 1:maxit){

    #E step
    #     r <- mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
    #                      U_mu0=sapply(U_mu, "["),
    #                      U_kappa0=sapply(U_kappa, "["),
    #                      U_nu0=sapply(U_nu, "["),
    #                      U_lambda0=U_Sigma,Log=TRUE)


    r <- mmNiWpdfC(Mu=sapply(mu_list, "["), Sigma=S_list,
                   U_Mu0=sapply(U_mu, "["),
                   U_Kappa0=sapply(U_kappa, "["),
                   U_Nu0=sapply(U_nu, "["),
                   U_Sigma0=U_Sigma,Log=TRUE)


    logexptrick_const <- apply(X=r, MARGIN=2, FUN=max)
    r_log_const <- apply(X=r, MARGIN=1, FUN=function(rv){rv-logexptrick_const})
    rw_const <- apply(X=exp(r_log_const), MARGIN=1, FUN=function(rv){rv*weights})
    r <- t(apply(X=rw_const, MARGIN=1, FUN=function(rv){rv/colSums(rw_const)}))


    #M step
    N_k <- rowSums(r)
    if(any(is.nan(N_k)) | any(round(N_k)==1) | any(!N_k) |any(rowSums((r==0)*1)==length(r[1,])-1)){
      # the above line is because sometimes for a given group, there is only one (or zero) loglikelihood superior to the likelihoods
      # computed for the other groups
      for(k in 1:K){
        #sampling the cluster parameters

        NiW <- rNiW(hyperG0, diagVar=FALSE)
        U_mu[[k]] <-NiW[["mu"]]
        U_Sigma[[k]] <-NiW[["S"]]
        U_nu[[k]]<-d+1
        U_kappa[[k]]<-0.01
      }
      maxit<-maxit+1
      next
    }
    weights  <- N_k/N
    #cat("weights:", weights, "\n")

    for(k in 1:K){

      rSinv_list <- mapply(S = S_list,
                           rik = as.list(r[k, ]),
                           FUN=function(S, rik){rik*solve(S)},
                           SIMPLIFY=FALSE)

      rSinv_sum <- Reduce('+', rSinv_list)


      U_mu[[k]]<-Reduce('+',lapply(1:N,
                                   function(i,m) {r[k,i]*mu_list[[i]]%*%m[[i]]},
                                   m=rSinv_list))%*% solve(rSinv_sum)


      mu_m <- lapply(mu_list, function(x){x - as.vector(U_mu[[k]])})

      U_kappa[[k]] <- (N_k[k]*d)/sum(mapply(x = mu_m, rSinv = rSinv_list,
                                            FUN=function(x,rSinv){
                                              t(x)%*%rSinv%*%x
                                            }, SIMPLIFY=TRUE))

      const_nu0_uniroot <-(sum(r[k,]*sapply(S_list, function(S){log(det(S))}))
                           + N_k[k]*log(det(rSinv_sum)))


      max_upper<-1e12
      U_nu[[k]] <- try(stats::uniroot(function(nu0){(N_k[k]*digamma_mv(x=nu0/2, p=d)
                                                     - N_k[k]*d*log(N_k[k]*nu0/2)
                                                     + const_nu0_uniroot
      )}, lower = d+1, upper=max_upper)$root, TRUE)


      if(inherits(U_nu[[k]], "try-error")){U_nu[[k]] <- d+1}
      U_nu[[k]]<-round(U_nu[[k]])  # warning issued in wishrnd: n isn't an integer anymore


      U_Sigma[[k]] <- N_k[k]*U_nu[[k]]*solve(rSinv_sum)
    }



    loglik[i+1] <- sum(log(colSums(weights%*%mmNiWpdf(mu=sapply(mu_list, "["), Sigma=S_list,
                                                      U_mu0=sapply(U_mu, "["),
                                                      U_kappa0=sapply(U_kappa, "["),
                                                      U_nu0=sapply(U_nu, "["),
                                                      U_lambda0=U_Sigma,Log=FALSE))))


    if (is.nan(loglik[i+1])){break}
    #cat("it ", i, ": loglik = ", loglik[i+1],"\n\n", sep="")
    if(abs(loglik[i+1]-loglik[i])<tol){break}
    #if (loglik[i+1]<loglik[i]){browser()}
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
              "U_mu" = U_mu,
              "U_kappa"=U_kappa,
              "U_nu" = U_nu,
              "U_lambda" = U_Sigma,
              "weights"=weights))

}
