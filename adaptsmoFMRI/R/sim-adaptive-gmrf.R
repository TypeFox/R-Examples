#' This function estimates the effects of a synthetic spatiotemporal data set resembling functional MR Images (fMRI), with the method of efficient Markov Chain Monte  Carlo (MCMC) simulation. The Metropolis Hastings (MH) algorithm is used for the non-approximate case and the Gibbs sampler for the approximate case.
#'
#' 
#'
#' @name sim.adaptiveGMRF
#' @aliases sim.adaptiveGMRF
#' @title Adaptive GMRF Model for Simulated Data
#' @usage sim.adaptiveGMRF(data, hrf, approximate = FALSE, K = 500, 
#' a = 1, b = 1, c = 1, d = 1, nu = 1, block = 1, burnin = 1, thin = 1)
#' @param data simulated fMRI-data, needs to be an array of dimension \code{(20 x 20 x T)}.
#' @param hrf haemodynamic response function, needs to be a vector of length \code{T}.
#' @param approximate logical, if \code{TRUE} then the approximate case is chosen. Default is \code{FALSE}.
#' @param K scalar, length of the MCMC path, hence iteration steps.
#' @param a scalar, shape hyperparameter of the inverse-gamma distribution of the variance parameter (\eqn{\sigma_i^2}).
#' @param b scalar, scale hyperparameter of the inverse gamma distribution of the variance parameter (\eqn{\sigma_i^2}).
#' @param c scalar, shape hyperparameter of the inverse gamma distribution of the precision parameter (\eqn{\tau}).
#' @param d scalar, scale hyperparameter of the inverse gamma distribution of the precision parameter (\eqn{\tau}).
#' @param nu scalar, shape and scale hyperparameter of the gamma distribution of the interaction weights (\eqn{w_{ij}}).
#' @param block scalar, when \code{approximate==TRUE} then a block of weights is updated at a time.
#' @param burnin scalar, defining the first iteration steps which should be omitted from MCMC path.
#' @param thin scalar, only every \code{thin} step of MCMC path is saved to output.
#' @author Max Hughes
#' @note This function is solely for one covariate.
#' @examples  
#' # non-transformed hr-function
#' T <- 210
#' seq.length <- T*3
#' index <- seq(3, T*3, by = 3)
#' hrf <- rep(c(-0.5, 0.5), each=30, times=ceiling(T/30*1.5))
#' hrf <- as.matrix(hrf[index])  
#' # get simulated data
#' data("sim_fmri")
#' data <- data_simfmri
#' # execute function
#' set.seed(111222)
#' K <- 2
#' a <- b <- c <- d <- nu <- 1
#' test.sim.adaptive <- sim.adaptiveGMRF(data, hrf, approximate=TRUE, K, 
#'                                       a, b, c, d, nu)


sim.adaptiveGMRF <- function(data, hrf, approximate=FALSE, K=500, 
                             a=1, b=1, c=1, d=1, nu=1, block=1, burnin=1, thin=1){

  ## load required libraries
  require("coda")
  require("mvtnorm")
  require("MCMCpack")
  require("Matrix")
  
  if(any(is.na(data)))
    stop("\nNAs in fMRI data.\n")
  
  if(dim(data)[1]!=20 || dim(data)[2]!=20 || (dim(data)[3]<60 && dim(data)[3]>900))
    stop("FMRI data needs to be an array of dimension (20x20x60<T<900).")

  if(any(is.na(hrf)))
    stop("\nNAs in hr function.\n")
  
  Z <- as.matrix(hrf)
  if(nrow(Z)!=dim(data)[3])
    stop("Haemodynamic response function needs to be of same length T (time) as
          simulated fMRI data.")
  
  if(a <= 0 || b <= 0 || c <= 0 || d <= 0)
    stop("Scale and shape parameters of the inverse gamma distributions need to be > 0.")
  
  dx <- dim(data)[1]
  dy <- dim(data)[2]
  I <- dx*dy
  T <- dim(data)[3]
  
  
  ## coerce y array with dim (nrow x ncol x T) to array
  y <- c()
  for (i in 1:dx){ 
    for (j in 1:dy){                
      y <- c(y, data[i,j,])  
    }
  }
 
  ## build K, as sparse matrix (adopted from cmr_space.R)
  # get coordinates
  coord <- c()
  for(i in 1:dx){
    for(j in 1:dy){
      coord <- cbind(coord, c(i, j))
    }
  }   #plot(as.matrix(coord)[1,], as.matrix(coord)[2,])

  nei <- c()
  for(i in 1:(I-1)){
    for(j in (i+1):I){
      if(sum((coord[,i]-coord[,j])^2)<2){
        nei <- cbind(nei, c(i, j))
      }
    }
  }   #plot(as.matrix(nei)[1,], as.matrix(nei)[2,])

  NEI <- dim(nei)[2]
  
  w <- rep(0.8, NEI)
  K.i <- c(1:I, nei[1,], nei[2,])
  K.j <- c(1:I, nei[2,], nei[1,])
  tauk.sq2K <- array(0, dim=c(length(K.i), NEI))

  for(i in 1:NEI){
    tauk.sq2K[nei[1, i], i] <- w[i]
    tauk.sq2K[nei[2, i], i] <- w[i]
    tauk.sq2K[I+i,i] <- - w[i]
    tauk.sq2K[I+NEI+i,i] <- - w[i]
  }   #image(as(tauk.sq2K, "sparseMatrix"))

  tauk.sq <- rep(1, NEI)
  K.sparse <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K%*%tauk.sq), dims=c(I,I))

  ## starting values
  beta <- numeric(I)
  w.new <- numeric(block)
  diff <- NEI - floor(NEI/block)*block
  sigma.sq <- rep(1, I)
  
  count <- 0
  acc.count <- 0
  
  ## Prepare mean e and precision Q for Rue/Held Algorithm
  tZZ <- t(Z)%*%Z
  sigtZZ <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ)
  sigtZy <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z))%*%y

  ## save output of MCMC
  beta.out <- w.out <- tauk.out <- sigma.out <- c()

  for(k in 1:K){  

    
    ## Step 1: Draw the blocks alpha_k from the Gaussian full    
    ##         conditionals k = 1,...,m.

    
    #             --  in simulated fmri data missing  --


    
    ## Step 2: Draw the blocks beta_k from the (multivariate) Gaussian full
    ##         conditionals
    
    # update Q
    Q <- sigtZZ + K.sparse
    # update e
    e <- sigtZy
    
    ## Rue/Held Algorithm, sampling beta ~ N(b, Q):
    # Step 1
    L <- chol(Q)  #as.matrix()
    # Step 2
    s <- solve(L,e)
    # Step 3
    mu.beta <- solve(t(L), s)
    # Step 4
    z <- rnorm(I)
    # Step 5
    v <- solve(t(L), z)
    # Step 6
    beta <- mu.beta+v
    

    ## Step 3: Draw the weights w_ij via MH steps or Gibbs-Sampling in the
    ##         approximate case
 
    if(approximate==TRUE){
      f.full <- nu/2
      for(i in 1:NEI){
        g.full <- nu/2+((beta[K.i[I+i]] - beta[K.j[I+i]])^2)/(2*tauk.sq[i])
        w[i] <- rgamma(1, shape=f.full, rate=g.full)
        # update tauk.sq2K
        tauk.sq2K[nei[1, i], i] <- w[i]
        tauk.sq2K[nei[2, i], i] <- w[i]
        tauk.sq2K[I+i,i] <- -w[i]
        tauk.sq2K[I+NEI+i,i] <- -w[i]
      }
    }
    else{
      K.old <- K.sparse
      f.full <- nu/2
      for(i in 1:NEI){
        # draw a block of weights at a time
        if(i%%block==0){
          count <- count + 1

          for(j in (1+i-block):i){
            g.full <- nu/2 + ((beta[K.i[I+j]] - beta[K.j[I+j]])^2)/(2*tauk.sq[j])
            w.new[j] <- rgamma(1, shape=f.full, rate=g.full)
            # update tauk.sq2K
            tauk.sq2K[nei[1,j],j] <- w.new[j]
            tauk.sq2K[nei[2,j],j] <- w.new[j]
            tauk.sq2K[I+j,j] <- -w.new[j]
            tauk.sq2K[I+NEI+j,j] <- -w.new[j]
          }
          # update K.sparse with w.new and call it K.new 
          K.new <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K%*%(1/tauk.sq)),
                                dims=c(I,I))

          # start the Cholesky decomposition in the row corresponding to the position
          # where a proposed new weight is located
          # delete last row and column
          eigen.new <- sum(log(diag(chol(K.new[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))
          eigen.old <- sum(log(diag(chol(K.old[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))
          
          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old==0||(eigen.old==Inf & eigen.new==Inf)){
            acc.rate <- 1
          }
          else{
            acc.rate <- sqrt(exp(eigen.new-eigen.old))
          }

          # accept
          if(runif(1)<acc.rate){
            acc.count <- acc.count + 1
            w[(1+i-block):i] <- w.new[(1+i-block):i]
          }
          else{
            w[(1+i-block):i] <- w[(1+i-block):i]
          }

          # update tauk.sq2K
          for(j in (1+i-block):i){
            tauk.sq2K[nei[1, j], j] <- w[j]
            tauk.sq2K[nei[2, j], j] <- w[j]
            tauk.sq2K[I+j,j] <- -w[j]
            tauk.sq2K[I+NEI+j,j] <- -w[j]
          }
        }
        # handle the weights which are left over
        if(i==NEI & NEI%%block!=0){
          count <- count + 1
          for(j in (1+i-diff):i){
            g.full <- nu/2 + ((beta[K.i[I+j]] - beta[K.j[I+j]])^2)/(2*tauk.sq[j])
            w.new[j] <- rgamma(1, shape=f.full, rate=g.full)
            # update tauk.sq2K
            tauk.sq2K[nei[1,j],j] <- w.new[j]
            tauk.sq2K[nei[2,j],j] <- w.new[j]
            tauk.sq2K[I+j,j] <- -w.new[j]
            tauk.sq2K[I+NEI+j,j] <- -w.new[j]
          }
          # update K.sparse with w.new and call it K.new 
          K.new <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K%*%(1/tauk.sq)),
                                dims=c(I,I))
          
          eigen.new <- sum(log(diag(chol(K.new[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))
 
          eigen.old <- sum(log(diag(chol(K.old[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))
     
          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old==0||(eigen.old==Inf & eigen.new==Inf)){
            acc.rate <- 1
          }
          else{
            acc.rate <- sqrt(exp(eigen.new-eigen.old))
          }

          # accept
          if(runif(1)<acc.rate){
            acc.count <- acc.count + 1
            w[(1+i-diff):i] <- w.new[(1+i-diff):i]
          }
          else{
            w[(1+i-diff):i] <- w[(1+i-diff):i]
          }      
          # update tauk.sq2K
          for(j in (1+i-diff):i){
            tauk.sq2K[nei[1,j],j] <- w[j]
            tauk.sq2K[nei[2,j],j] <- w[j]
            tauk.sq2K[I+j,j] <- -w[j]
            tauk.sq2K[I+NEI+j,j] <- -w[j]
          }
        }
      } 
    }        

      
    ### Step 4: Draw the variance parameters sigma^2 and the hyperparameters tau^2 from
    ###         their corresponding inverse gamma full conditionals.

    ## sigma^2
    a.full <- a + 0.5*T
    for(i in 1:I){
      b.full <- b + 0.5*sum((y[1:T+(i-1)*T] - (Z%*%beta[i]))^2)  
      sigma.sq[i] <- rinvgamma(1, a.full, b.full)
    }

    # update sigtZZ
    sigtZZ <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ)
    # update sigtZy
    sigtZy <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z))%*%y
  
    # tau^2
    c.full <- c + 0.5*(I-1)
    d.full <- d + 0.5*t(beta)%*%K.sparse%*%beta
    tauk.sq <- rinvgamma(NEI, c.full, d.full[1,1])
    tauk.sq <- rep(median(tauk.sq), NEI)

    # update K.sparse
    K.sparse <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K%*%(1/tauk.sq)),
                             dims=c(I, I))
    
    print(c(k, median(beta), median(w), median(sigma.sq), median(tauk.sq)))
    
    ## start saving MCMC output after burnin
    if(k >= burnin){
      ## save only every thin iteration of MCMC output
      if(k%%thin==0){
        beta.out <- cbind(beta.out, as.vector(beta))
        w.out <- cbind(w.out, w)
        tauk.out <- c(tauk.out, median(tauk.sq))
        sigma.out <- cbind(sigma.out, sigma.sq)
      }
    }
  }

  result <- list("dx"=dx, "dy"=dy, "I"=I, "iter"=K, "coord"=coord, "nei"=nei,
                 "NEI"=NEI, "count"=count, "acc.count"=acc.count, 
                 "beta.out"=beta.out, "w.out"=w.out,
                 "sigma.out"=sigma.out, "tauk.out"=tauk.out)
  class(result) <- "adaptsmoFMRI"
  return(result)
}

