#' This function estimates the effects of a synthetic spatiotemporal data set resembling functional MR Images (fMRI), with the method of efficient Markov Chain Monte  Carlo (MCMC) simulation. The Metropolis Hastings (MH) algorithm is used for the non-approximate case and the Gibbs sampler for the approximate case.
#'
#' 
#'
#' @name sim.adaptiveGMRF2COVAR
#' @aliases sim.adaptiveGMRF2COVAR
#' @title Adaptive GMRF Model for Simulated Data
#' @usage sim.adaptiveGMRF2COVAR(data, hrf, approximate = FALSE, K = 500,
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
#' @note This function is solely for two covariates.
#' @examples
#' # See example function for simulated data (one covariate).        


sim.adaptiveGMRF2COVAR <- function(data, hrf, approximate=FALSE, K=500,
                                   a=1, b=1, c=1, d=1, nu=1, block=1, burnin=1, thin=1){

  ## load required libraries
  require("coda")
  require("mvtnorm")
  require("MCMCpack")
  require("Matrix")
  require("parallel")
  
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

  ## For model with two covariates
  p <- dim(Z)[2]
  if(p!=2)
    stop("Haemodynamic response function needs to be a matrix with column dimension of 2.")
    
  if(a <= 0 || b <= 0 || c <= 0 || d <= 0)
    stop("Scale and shape parameters of the inverse gamma distributions need to be > 0.")

  dx <- dim(data)[1]
  dy <- dim(data)[2]
  I <- dx*dy
  T <- dim(data)[3]
  Z.Var1 <- as.matrix(Z[,1])
  Z.Var2 <- as.matrix(Z[,2])
  
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

  K.i <- c(1:I, nei[1,], nei[2,])
  K.j <- c(1:I, nei[2,], nei[1,])
  
  ## K matrix for first covariate
  w.Var1 <- rep(0.8, NEI)
  tauk.sq2K.Var1 <- as(array(0, dim=c(length(K.i), NEI)), "sparseMatrix")

  for(i in 1:NEI){
    tauk.sq2K.Var1[nei[1, i], i] <- w.Var1[i]
    tauk.sq2K.Var1[nei[2, i], i] <- w.Var1[i]
    tauk.sq2K.Var1[I+i,i] <- -w.Var1[i]
    tauk.sq2K.Var1[I+NEI+i,i] <- -w.Var1[i]
  }   #image(as(tauk.sq2K.Var1, "sparseMatrix"))


  ## precision paramter for first covariate
  tauk.sq.Var1 <- rep(1, NEI)
  ## tauk.sq.Var1 togehter with tauk.sq.Var1
  K.sparse.Var1 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var1%*%tauk.sq.Var1), dims=c(I,I))

  ## K matrix for second covariate
  w.Var2 <- rep(0.8, NEI)
  tauk.sq2K.Var2 <- as(array(0, dim=c(length(K.i), NEI)), "sparseMatrix")

  for(i in 1:NEI){
    tauk.sq2K.Var2[nei[1, i], i] <- w.Var2[i]
    tauk.sq2K.Var2[nei[2, i], i] <- w.Var2[i]
    tauk.sq2K.Var2[I+i,i] <- -w.Var2[i]
    tauk.sq2K.Var2[I+NEI+i,i] <- -w.Var2[i]
  }   #image(as(tauk.sq2K.Var2, "sparseMatrix"))

  ## precision paramter for first covariate
  tauk.sq.Var2 <- rep(1, NEI)
  ## tauk.sq.Var2 togehter with tauk.sq.Var2
  K.sparse.Var2 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var2%*%tauk.sq.Var2), dims=c(I,I))
  
  ## starting values
  # Beta for first covariate
  beta.Var1 <- numeric(I)
  # Beta for second covariate
  beta.Var2 <- numeric(I)
  # proposed new weights for first covariate
  w.new.Var1 <- numeric(block)
  # proposed new weights for second covariate
  w.new.Var2 <- numeric(block)
  diff <- NEI - floor(NEI/block)*block
  sigma.sq <- rep(1, I)

  ## count for Variable 1
  count.Var1 <- 0
  acc.count.Var1 <- 0
  
  ## count for Variable 2
  count.Var2 <- 0
  acc.count.Var2 <- 0

  ## Prepare mean e.Var1 and precision Q.Var1 for Rue/Held Algorithm
  # First variable
  tZZ.Var1 <- t(Z.Var1)%*%Z.Var1
  sigtZZ.Var1 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ.Var1)
  sigtZy.Var1 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z.Var1))%*%y
                                                      #here beta.Var2 = 0
  #Second variable
  ## Prepare mean e.Var2 and precision Q.Var2 for Rue/Held Algorithm
  tZZ.Var2 <- t(Z.Var2)%*%Z.Var2
  sigtZZ.Var2 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ.Var2)
  sigtZy.Var2 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z.Var2))%*%y
                                                      #here beta.Var1 = 0

  ## save output of MCMC
  beta.Var1.out <- beta.Var2.out <- w.Var1.out <- w.Var2.out <- tauk.Var1.out <- tauk.Var2.out <- sigma.out <- c()

  for(k in 1:K){


    ## Step 1: Draw the blocks alpha_k from the Gaussian full
    ##         conditionals k = 1,...,m.


    #             --  in simulated fmri data missing  --



    ## Step 2: Draw the blocks beta_k from the (multivariate) Gaussian full
    ##         conditionals

    ## Vaiable 1
    # update Q
    Q.Var1 <- sigtZZ.Var1 + K.sparse.Var1
    # update e
    e.Var1 <- sigtZy.Var1

    ## Rue/Held Algorithm, sampling beta ~ N(b, Q):
    # Step 1
    L.Var1 <- chol(Q.Var1)  #as.matrix()
    # Step 2
    s.Var1 <- solve(L.Var1,e.Var1)
    # Step 3
    mu.beta.Var1 <- solve(t(L.Var1), s.Var1)
    # Step 4
    z.Var1 <- rnorm(I)
    # Step 5
    v.Var1 <- solve(t(L.Var1), z.Var1)
    # Step 6
    beta.Var1 <- mu.beta.Var1+v.Var1

    ## Vaiable 2
    # update Q
    Q.Var2 <- sigtZZ.Var2 + K.sparse.Var2
    # update e
    e.Var2 <- sigtZy.Var2

    ## Rue/Held Algorithm, sampling beta ~ N(b, Q):
    # Step 1
    L.Var2 <- chol(Q.Var2)  #as.matrix()
    # Step 2
    s.Var2 <- solve(L.Var2,e.Var2)
    # Step 3
    mu.beta.Var2 <- solve(t(L.Var2), s.Var2)
    # Step 4
    z.Var2 <- rnorm(I)
    # Step 5
    v.Var2 <- solve(t(L.Var2), z.Var2)
    # Step 6
    beta.Var2 <- mu.beta.Var2+v.Var2


    ## Step 3: Draw the weights w_ij via MH steps or Gibbs-Sampling in the
    ##         approximate case

    if(approximate==TRUE){
      ## Update weights for first covariate
      f.full.Var1 <- nu/2
      for(i in 1:NEI){
        g.full.Var1 <- nu/2+((beta.Var1[K.i[I+i]] - beta.Var1[K.j[I+i]])^2)/(2*tauk.sq.Var1[i])
        w.Var1[i] <- rgamma(1, shape=f.full.Var1, rate=g.full.Var1)
        # update tauk.sq2K
        tauk.sq2K.Var1[nei[1, i], i] <- w.Var1[i]
        tauk.sq2K.Var1[nei[2, i], i] <- w.Var1[i]
        tauk.sq2K.Var1[I+i,i] <- -w.Var1[i]
        tauk.sq2K.Var1[I+NEI+i,i] <- -w.Var1[i]
      }
      ## Update weights for second covariate
      f.full.Var2 <- nu/2
      for(i in 1:NEI){
        g.full.Var2 <- nu/2+((beta.Var2[K.i[I+i]] - beta.Var2[K.j[I+i]])^2)/(2*tauk.sq.Var2[i])
        w.Var2[i] <- rgamma(1, shape=f.full.Var2, rate=g.full.Var2)
        # update tauk.sq2K
        tauk.sq2K.Var2[nei[1, i], i] <- w.Var2[i]
        tauk.sq2K.Var2[nei[2, i], i] <- w.Var2[i]
        tauk.sq2K.Var2[I+i,i] <- -w.Var2[i]
        tauk.sq2K.Var2[I+NEI+i,i] <- -w.Var2[i]
      }
    }
    else{
      ## Update weights for first covariate
      drawWeights.Var1 <- function(...){
      K.old.Var1 <- K.sparse.Var1
      f.full.Var1 <- nu/2
      for(i in 1:NEI){
        # draw a block of weights at a time
        if(i%%block==0){
          count.Var1 <- count.Var1 + 1
          for(j in (1+i-block):i){
            g.full.Var1 <- nu/2 + ((beta.Var1[K.i[I+j]] - beta.Var1[K.j[I+j]])^2)/(2*tauk.sq.Var1[j])
            w.new.Var1[j] <- rgamma(1, shape=f.full.Var1, rate=g.full.Var1)
            # update tauk.sq2K
            tauk.sq2K.Var1[nei[1,j],j] <- w.new.Var1[j]
            tauk.sq2K.Var1[nei[2,j],j] <- w.new.Var1[j]
            tauk.sq2K.Var1[I+j,j] <- -w.new.Var1[j]
            tauk.sq2K.Var1[I+NEI+j,j] <- -w.new.Var1[j]
          }
          # update K.sparse with w.new and call it K.new
          K.new.Var1 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var1%*%(1/tauk.sq.Var1)),
                                dims=c(I,I))
          # start the Cholesky decomposition in the row corresponding to the position
          # where a proposed new weight is located
          # delete last row and column
          eigen.new.Var1 <- sum(log(diag(chol(K.new.Var1[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))
          eigen.old.Var1 <- sum(log(diag(chol(K.old.Var1[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))

          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old.Var1==0||(eigen.old.Var1==Inf & eigen.new.Var1==Inf)){
            acc.rate.Var1 <- 1
          }
          else{
            acc.rate.Var1 <- sqrt(exp(eigen.new.Var1-eigen.old.Var1))
          }

          # accept
          if(runif(1)<acc.rate.Var1){
            acc.count.Var1 <- acc.count.Var1 + 1
            w.Var1[(1+i-block):i] <- w.new.Var1[(1+i-block):i]
          }
          else{
            w.Var1[(1+i-block):i] <- w.Var1[(1+i-block):i]
          }

          # update tauk.sq2K
          for(j in (1+i-block):i){
            tauk.sq2K.Var1[nei[1, j], j] <- w.Var1[j]
            tauk.sq2K.Var1[nei[2, j], j] <- w.Var1[j]
            tauk.sq2K.Var1[I+j,j] <- -w.Var1[j]
            tauk.sq2K.Var1[I+NEI+j,j] <- -w.Var1[j]
          }
        }
        # handle the weights which are left over
        if(i==NEI & NEI%%block!=0){
          count.Var1 <- count.Var1 + 1
          for(j in (1+i-diff):i){
            g.full.Var1 <- nu/2 + ((beta.Var1[K.i[I+j]] - beta.Var1[K.j[I+j]])^2)/(2*tauk.sq.Var1[j])
            w.new.Var1[j] <- rgamma(1, shape=f.full.Var1, rate=g.full.Var1)
            # update tauk.sq2K
            tauk.sq2K.Var1[nei[1,j],j] <- w.new.Var1[j]
            tauk.sq2K.Var1[nei[2,j],j] <- w.new.Var1[j]
            tauk.sq2K.Var1[I+j,j] <- -w.new.Var1[j]
            tauk.sq2K.Var1[I+NEI+j,j] <- -w.new.Var1[j]
          }
          # update K.sparse with w.new and call it K.new
          K.new.Var1 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var1%*%(1/tauk.sq.Var1)),
                                dims=c(I,I))

          eigen.new.Var1 <- sum(log(diag(chol(K.new.Var1[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))

          eigen.old.Var1 <- sum(log(diag(chol(K.old.Var1[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))

          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old.Var1==0||(eigen.old.Var1==Inf & eigen.new.Var1==Inf)){
            acc.rate.Var1 <- 1
          }
          else{
            acc.rate.Var1 <- sqrt(exp(eigen.new.Var1-eigen.old.Var1))
          }

          # accept
          if(runif(1)<acc.rate.Var1){
            acc.count.Var1 <- acc.count.Var1 + 1
            w.Var1[(1+i-diff):i] <- w.new.Var1[(1+i-diff):i]
          }
          else{
            w.Var1[(1+i-diff):i] <- w.Var1[(1+i-diff):i]
          }
          # update tauk.sq2K
          for(j in (1+i-diff):i){
            tauk.sq2K.Var1[nei[1,j],j] <- w.Var1[j]
            tauk.sq2K.Var1[nei[2,j],j] <- w.Var1[j]
            tauk.sq2K.Var1[I+j,j] <- -w.Var1[j]
            tauk.sq2K.Var1[I+NEI+j,j] <- -w.Var1[j]
          }
        }
      }
      
      return(list("w.Var1"=w.Var1, "tauk.sq2K.Var1"=tauk.sq2K.Var1,
                  "count.Var1"=count.Var1, "acc.count.Var1"=acc.count.Var1))
    }
      ## Update weights for second covariate
      drawWeights.Var2 <- function(...){
      K.old.Var2 <- K.sparse.Var2
      f.full.Var2 <- nu/2
      for(i in 1:NEI){
        # draw a block of weights at a time
        if(i%%block==0){
          count.Var2 <- count.Var2 + 1
          for(j in (1+i-block):i){
            g.full.Var2 <- nu/2 + ((beta.Var2[K.i[I+j]] - beta.Var2[K.j[I+j]])^2)/(2*tauk.sq.Var2[j])
            w.new.Var2[j] <- rgamma(1, shape=f.full.Var2, rate=g.full.Var2)
            # update tauk.sq2K
            tauk.sq2K.Var2[nei[1,j],j] <- w.new.Var2[j]
            tauk.sq2K.Var2[nei[2,j],j] <- w.new.Var2[j]
            tauk.sq2K.Var2[I+j,j] <- -w.new.Var2[j]
            tauk.sq2K.Var2[I+NEI+j,j] <- -w.new.Var2[j]
          }
          # update K.sparse with w.new and call it K.new
          K.new.Var2 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var2%*%(1/tauk.sq.Var2)),
                                dims=c(I,I))

          # start the Cholesky decomposition in the row corresponding to the position
          # where a proposed new weight is located
          # delete last row and column
          eigen.new.Var2 <- sum(log(diag(chol(K.new.Var2[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))
          eigen.old.Var2 <- sum(log(diag(chol(K.old.Var2[K.i[I+1-block+i]:(I-1),
                                   K.i[I+1-block+i]:(I-1)]))))

          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old.Var2==0||(eigen.old.Var2==Inf & eigen.new.Var2==Inf)){
            acc.rate.Var2 <- 1
          }
          else{
            acc.rate.Var2 <- sqrt(exp(eigen.new.Var2-eigen.old.Var2))
          }

          # accept
          if(runif(1)<acc.rate.Var2){
            acc.count.Var2 <- acc.count.Var2 + 1
            w.Var2[(1+i-block):i] <- w.new.Var2[(1+i-block):i]
          }
          else{
            w.Var2[(1+i-block):i] <- w.Var2[(1+i-block):i]
          }

          # update tauk.sq2K
          for(j in (1+i-block):i){
            tauk.sq2K.Var2[nei[1, j], j] <- w.Var2[j]
            tauk.sq2K.Var2[nei[2, j], j] <- w.Var2[j]
            tauk.sq2K.Var2[I+j,j] <- -w.Var2[j]
            tauk.sq2K.Var2[I+NEI+j,j] <- -w.Var2[j]
          }
        }
        # handle the weights which are left over
        if(i==NEI & NEI%%block!=0){
          count.Var2 <- count.Var2 + 1
          for(j in (1+i-diff):i){
            g.full.Var2 <- nu/2 + ((beta.Var2[K.i[I+j]] - beta.Var2[K.j[I+j]])^2)/(2*tauk.sq.Var2[j])
            w.new.Var2[j] <- rgamma(1, shape=f.full.Var2, rate=g.full.Var2)
            # update tauk.sq2K
            tauk.sq2K.Var2[nei[1,j],j] <- w.new.Var2[j]
            tauk.sq2K.Var2[nei[2,j],j] <- w.new.Var2[j]
            tauk.sq2K.Var2[I+j,j] <- -w.new.Var2[j]
            tauk.sq2K.Var2[I+NEI+j,j] <- -w.new.Var2[j]
          }
          # update K.sparse with w.new and call it K.new
          K.new.Var2 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var2%*%(1/tauk.sq.Var2)),
                                dims=c(I,I))

          eigen.new.Var2 <- sum(log(diag(chol(K.new.Var2[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))

          eigen.old.Var2 <- sum(log(diag(chol(K.old.Var2[K.i[I+1-block+i]:(I-1),
                                    K.i[I+1-block+i]:(I-1)]))))

          # avoid deviding by 0 and avoid having infinity in numerator and denominator!
          # (better set acc.rate to zero in infinity case?)
          if(eigen.old.Var2==0||(eigen.old.Var2==Inf & eigen.new.Var2==Inf)){
            acc.rate.Var2 <- 1
          }
          else{
            acc.rate.Var2 <- sqrt(exp(eigen.new.Var2-eigen.old.Var2))
          }

          # accept
          if(runif(1)<acc.rate.Var2){
            acc.count.Var2 <- acc.count.Var2 + 1
            w.Var2[(1+i-diff):i] <- w.new.Var2[(1+i-diff):i]
          }
          else{
            w.Var2[(1+i-diff):i] <- w.Var2[(1+i-diff):i]
          }
          # update tauk.sq2K
          for(j in (1+i-diff):i){
            tauk.sq2K.Var2[nei[1,j],j] <- w.Var2[j]
            tauk.sq2K.Var2[nei[2,j],j] <- w.Var2[j]
            tauk.sq2K.Var2[I+j,j] <- -w.Var2[j]
            tauk.sq2K.Var2[I+NEI+j,j] <- -w.Var2[j]
          }
        }
      }

      return(list("w.Var2"=w.Var2, "tauk.sq2K.Var2"=tauk.sq2K.Var2,
                  "count.Var2"=count.Var2, "acc.count.Var2"=acc.count.Var2))
    }
      weights.Var1 <- mcparallel(drawWeights.Var1(K.sparse.Var1, nu, nei, NEI, I, block, count.Var1, beta.Var1, tauk.sq.Var1, tauk.sq2K.Var1, K.i, K.j), name="Var1")
      weights.Var2 <- mcparallel(drawWeights.Var2(K.sparse.Var2, nu, nei, NEI, I, block, count.Var2, beta.Var2, tauk.sq.Var2, tauk.sq2K.Var2, K.i, K.j), name="Var2")
      results <- mccollect(list("Var1"=weights.Var1, "Var2"=weights.Var2))

      w.Var1 <- results$Var1$w.Var1
      tauk.sq2K.Var1 <- results$Var1$tauk.sq2K.Var1
      count.Var1 <- results$Var1$count.Var1
      acc.count.Var1 <- results$Var1$acc.count.Var1
      
      w.Var2 <- results$Var2$w.Var2
      tauk.sq2K.Var2 <- results$Var2$tauk.sq2K.Var2
      count.Var2 <- results$Var2$count.Var2
      acc.count.Var2 <- results$Var2$acc.count.Var2
    }
      
    ### Step 4: Draw the variance parameters sigma^2 and the hyperparameters tau^2 from
    ###         their corresponding inverse gamma full conditionals.

    ## sigma^2
    a.full <- a + 0.5*T
    for(i in 1:I){
      b.full <- b + 0.5*sum((y[1:T+(i-1)*T] -
                             (Z.Var1%*%beta.Var1[i] +
                              Z.Var2%*%beta.Var2[i]))^2)
      sigma.sq[i] <- rinvgamma(1, a.full, b.full)
    }

    # update sigtZZ.Var1
    sigtZZ.Var1 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ.Var1)
    # update sigtZy.Var1
    sigtZy.Var1 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z.Var1))%*%(y -
                                  kronecker(as(diag(1, I), "sparseMatrix"), Z.Var2)%*%beta.Var2)

    # update sigtZZ.Var2
    sigtZZ.Var2 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), tZZ.Var2)
    # update sigtZy.Var2
    sigtZy.Var2 <- kronecker(as(diag(1/sigma.sq, I), "sparseMatrix"), t(Z.Var2))%*%(y -
                                  kronecker(as(diag(1, I), "sparseMatrix"), Z.Var1)%*%beta.Var1)

    # tau^2 for Variable 1
    c.full.Var1 <- c+0.5*(I-1)
    d.full.Var1 <- d+0.5*t(beta.Var1)%*%K.sparse.Var1%*%beta.Var1
    tauk.sq.Var1 <- rinvgamma(NEI, c.full.Var1, d.full.Var1[1,1])
    tauk.sq.Var1 <- rep(median(tauk.sq.Var1), NEI)

    # update K.sparse
    K.sparse.Var1 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var1%*%(1/tauk.sq.Var1)),
                             dims=c(I,I))

    # tau^2 for Variable 2
    c.full.Var2 <- c+0.5*(I-1)
    d.full.Var2 <- d+0.5*t(beta.Var2)%*%K.sparse.Var2%*%beta.Var2
    tauk.sq.Var2 <- rinvgamma(NEI, c.full.Var2, d.full.Var2[1,1])
    tauk.sq.Var2 <- rep(median(tauk.sq.Var2), NEI)

    # update K.sparse
    K.sparse.Var2 <- sparseMatrix(K.i, K.j, x=as.vector(tauk.sq2K.Var2%*%(1/tauk.sq.Var2)),
                             dims=c(I,I))

    print(c(k, median(beta.Var1), median(beta.Var2), median(w.Var1), median(w.Var2), median(sigma.sq),
            median(tauk.sq.Var1), median(tauk.sq.Var2)))

    ## start saving MCMC output after burnin
    if(k >= burnin){
      ## save only every thin iteration of MCMC output
      if(k%%thin==0){
        beta.Var1.out <- cbind(beta.Var1.out, as.vector(beta.Var1))
        beta.Var2.out <- cbind(beta.Var2.out, as.vector(beta.Var2))
        w.Var1.out <- cbind(w.Var1.out, w.Var1)
        w.Var2.out <- cbind(w.Var2.out, w.Var2)
        tauk.Var1.out <- c(tauk.Var1.out, median(tauk.sq.Var1))
        tauk.Var2.out <- c(tauk.Var2.out, median(tauk.sq.Var2))
        sigma.out <- cbind(sigma.out, sigma.sq)
      }
    }
  }

  result <- list("dx"=dx, "dy"=dy, "I"=I, "iter"=K, "coord"=coord, "nei"=nei,
                 "NEI"=NEI, "count.Var1"=count.Var1, "acc.count.Var1"=acc.count.Var1,
                 "count.Var2"=count.Var2, "acc.count.Var2"=acc.count.Var2,
                 "beta.Var1.out"=beta.Var1.out, "beta.Var2.out"=beta.Var2.out,
                 "w.Var1.out"=w.Var1.out, "w.Var2.out"=w.Var2.out,
                 "sigma.out"=sigma.out,
                 "tauk.Var1.out"=tauk.Var1.out, "tauk.Var2.out"=tauk.Var2.out)
                 
  class(result) <- "adaptsmoFMRI"
  return(result)
}

