#
# Tools for analysing the models learned with the
# functions in CCAGFA.R
#
# Contains functions:
#  CCAtrim()   : Cleans up the CCA/BIBFA model by dropping inactive components
#  GFAtrim()   : The same for GFA
#  CCAcorr()   : Computes the correlations of the components
#  CCApred()   : Make predictions for new data for CCA/BIBFA
#  GFApred()   : Make predictions for new data for GFA
#  CCAsample() : Generate data from the model for CCA/BIBFA
#  GFApred()   : Make predictions for new data for GFA
#

CCAtrim <- function(model,threshold=1e-3) {
  #
  # A wrapper for GFAtrim(), to be used for CCA/BIBFA models.
  #

  if(length(model$D)!=2) {
    print("Not a CCA/BIBFA model. Use GFAtrim() instead.")
  } else {
    return(GFAtrim(model,threshold))
  }
}
  
GFAtrim <- function(model,threshold=1e-3) {
  #
  # Clean up the model by removing the components that were
  # pushed to zero.
  #
  # Input:
  #  model    : A model trained with GFA()
  #  threshold: Components that explain less than threshold
  #             of the total variation are dropped
  #
  # Output:
  #  model    : The trimmed model (see GFA() for the contents)
  #

  M <- length(model$W)
  N <- dim(model$Z)[1]

  active <- matrix(1,M,model$K) 
  for(m in 1:M) {
    # The relative contribution to describing the total data
    # variation is D[m]/alpha[m,k] / (datavar[m] - D[m]/tau[m])
    residual <- model$datavar[m] - model$D[m]/model$tau[m]
    if(residual < 1e-2*model$datavar[m]) {  # Noise models almost everything
      print(paste("Warning, trimming a model for which the noise explains already",
            format(model$D[m]/model$tau[m]/model$datavar[m],digits=4),"percent of the total variation in data set",m,". Might result for zero components for that data set."))
      residual <- 1e-2*model$datavar[m]
    }

    active[m,which(model$D[m]/model$alpha[m,]<threshold*residual)] <- 0
  }

  keep <- which(colSums(active)>0)
  model$Z <- model$Z[,keep,drop=F]
  model$covZ <- model$covZ[keep,keep,drop=F]
  
  if(!is.null(dim(model$covW))) { #covW list not stored due to low.mem=TRUE
    model$covW <- list()
    for(m in 1:M)
      model$covW[[m]] <- (model$WW[[m]] - crossprod(model$W[[m]]))/model$D[m]
  }

  for(m in 1:M) {
    model$W[[m]] <- model$W[[m]][,keep,drop=F]
    model$covW[[m]] <- model$covW[[m]][keep,keep,drop=F]
    model$WW[[m]] <- model$WW[[m]][keep,keep,drop=F]
  }
  model$alpha <- model$alpha[,keep,drop=F]

  active <- active[,keep,drop=F]
  model$K <- length(keep)

  model$ZZ <- crossprod(model$Z) + N*model$covZ
  id <- rep(1,model$K)   # Vector of ones for fast matrix calculations
  for(m in 1:M) {
    for(k in 1:model$K) {
      if(active[m,k]==0) {
        model$W[[m]][,k] <- 0
      }
    }
    tmp <- 1/sqrt(model$alpha[m,])
    model$covW[[m]] <- 1/model$tau[m] * outer(tmp,tmp) *
      chol2inv(chol(outer(tmp,tmp)*model$ZZ + diag(1/model$tau[m],model$K)))

    # An alternative way; could be tried in case of 
    #   issues with numerical stability
    #eS <- eigen( outer( tmp, id )*model$ZZ*outer(id,tmp) + diag(1/model$tau[m],model$K) , symmetric=TRUE)
    #model$covW[[m]] <- 1/tau[m] * outer( tmp, id ) * tcrossprod( eS$vectors*outer( id, 1/eS$values), eS$vectors ) * outer(id,tmp)

    model$WW[[m]] <- crossprod(model$W[[m]]) + model$covW[[m]]*model$D[m]
  }

  model$active <- active
  model$trimmed <- TRUE

  return(model)
}
  
CCAcorr <- function(Y,model,threshold=1e-3) {
  #
  # Function for estimating correlations as in the regular CCA.
  # Note that the CCA model does not guarantee orthogonality
  # of the individual components in the same way as CCA, so
  # the results might not be directly interpretable in this
  # way.
  #
  # Inputs:
  #  Y         : Test data as in CCA()
  #  model     : A model learned with CCA()
  #  threshold : Components explaining at least threshold of
  #              the total variance in both data sets are shared
  #
  # Output:
  # A list containing:
  #  r         : Correlation coefficients for each latent dimension
  #  active    : An indicator saying which of the components are
  #              shared; the correlations for the non-shared components
  #              are not interesting
  #

  if(length(model$D)!=2) {
    print("Correlation only makes sense for CCA/BIBFA models.")
    return
  }
  
  # Find the latent variables conditional on each data
  pred1 <- CCApred(c(1,0),Y,model)
  pred2 <- CCApred(c(0,1),Y,model)

  # Find the correlation for each component
  r <- vector(length=model$K)
  for(k in 1:model$K) {
    r[k] <- cor(pred1$Z[,k],pred2$Z[,k])
  }
  active <- matrix(1,length(Y),model$K)
  for(m in 1:length(Y)) {
    residual <- model$datavar[m] - model$D[m]/model$tau[m]
    active[m,which(model$D[m]/model$alpha[m,] < threshold*residual)] <- 0
  }
  active <- apply(active,2,min)

  return(list(r=r,active=active))
}

CCApred <- function(pred,Y,model,sample=FALSE,nSample=100) {
  #
  # A wrapper for GFApred, to be used for CCA/BIBFA models.
  #

  if(length(model$D)!=2) {
    print("Not a CCA/BIBFA model. Use GFApred() instead.")
  } else {
    return(GFApred(pred,Y,model,sample,nSample))
  }
}

GFApred <- function(pred,Y,model,sample=FALSE,nSample=100) {
  #
  # Function for making predictions with the model. Gives the
  # mean prediction and the mean and covariance of the latent
  # variables. The predictive distribution itself does not have
  # a closed-form expression, so the function also allows drawing
  # samples from it.
  #
  # Inputs:
  #   pred:  Binary vector of length 2, indicating which of the
  #          two data sets have been observed. (1,0) indicates
  #          we observe the first data set and want to predict
  #          the values for the latter, and (0,1) does the opposite.
  #          Using (1,1) allows computing the latent variables
  #          for new test samples where both views are observed.
  #   Y   :  The test data as a list of length 2, given in the
  #          same format as for the function GFA(). The data
  #          matrix for the missing views can be anything, e.g.
  #          zeros, but it needs to exist
  #   model: A model learned from training data using GFA()
  #   sample: Should we sample observations from the full predictive
  #           distribution?
  #   nSample: How many samples to draw if sample==TRUE
  #
  #
  # Outputs:
  # A list containing:
  #   Y    : The mean predictions as list. Observed data sets are retained
  #          as they were.
  #   Z    : Mean latent variables of the test samples, given the observed
  #          data; N times K matrix
  #   covZ : Covariance of the latent variables; K times K matrix
  #   sam  : Samples drawn from the predictive distribution, only
  #          returned if sample==TRUE. A list of Z, W and Y.
  #          Z is nSample times N times K matrix of the samples values.
  #          W and Y are M-element lists where only the predicted
  #          views are included (to avoid storing nSample identical
  #          copies of the observed data), each being a multidimensional
  #          array of nSample times the size of W and Y, respectively.
  #
  
  tr <- which(pred==1)  # The observed data sets
  pr <- which(pred==0)  # The data sets that need to be predicted

  N <- nrow(Y[[tr[1]]])
  M <- length(model$D)
  
  if(!is.null(dim(model$covW))) { #covW list not stored due to low.mem=TRUE
    model$covW <- list()
    for(m in 1:M)
      model$covW[[m]] <- (model$WW[[m]] - crossprod(model$W[[m]]))/model$D[m]
  }
  
  # Estimate the covariance of the latent variables
  covZ <- diag(1,model$K)
  for(m in tr) {
    covZ <- covZ + model$tau[m]*model$WW[[m]]
  }

  # Estimate the latent variables
  eS <- eigen( covZ ,symmetric=TRUE)
  covZ <- tcrossprod( eS$vectors*outer(rep(1,model$K),1/eS$values), eS$vectors )
  Z <- matrix(0,N,model$K)
  for(m in tr) {
    Z <- Z + Y[[m]]%*%model$W[[m]]*model$tau[m]
  }
  Z <- Z%*%covZ

  # Add a tiny amount of noise on top of the latent variables,
  # to supress possible artificial structure in components that
  # have effectively been turned off
  Z <- Z + model$addednoise*matrix(rnorm(N*model$K,0,1),N,model$K) %*% chol(covZ)

  # The prediction
  # NOTE: The ICML'11 paper has a typo in the prediction formula
  # on page 5. The mean prediction should have W_2^T instead of W_2.
  for(m in pr) {
    Y[[m]] <- tcrossprod(Z,model$W[[m]])
  }

  # Sample from the predictive distribution
  # Note that this code is fairly slow fow large nSample
  if(sample) {
    sam <- list()
    sam$Z <- array(0,c(nSample,N,model$K))
    sam$Y <- vector("list",length=M)
    sam$W <- vector("list",length=M)
    cholW <- vector("list",length=M)
    for(m in pr) {
      cholW[[m]] <- chol(model$covW[[m]])
      sam$W[[m]] <- array(0,c(nSample,model$D[m],model$K))
      sam$Y[[m]] <- array(0,c(nSample,N,model$D[m]))
    }

    cholZ <- chol(covZ)
    for(i in 1:nSample) {
      Ztemp <- Z + matrix(rnorm(N*model$K,0,1),N,model$K) %*% cholZ
      sam$Z[i,,] <- Ztemp
      for(m in pr) {
        Wtemp <- model$W[[m]] + matrix(rnorm(model$K*model$D[[m]],0,1),model$D[m],model$K) %*% cholW[[m]]
        sam$W[[m]][i,,] <- Wtemp
        sam$Y[[m]][i,,] <- tcrossprod(Ztemp,Wtemp) + matrix(rnorm(N*model$D[m],0,1/sqrt(model$tau[m])),N,model$D[m])
      }
    }
  }
  
  if(sample)
    return(list(Y=Y,Z=Z,covZ=covZ,sam=sam))
  else
    return(list(Y=Y,Z=Z,covZ=covZ))
}

CCAsample <- function(model,N) {
  #
  # A wrapper for GFApred, to be used for CCA/BIBFA models.
  #

  if(length(model$D)!=2) {
    print("Not a CCA/BIBFA model. Use GFAsample() instead.")
  } else {
    return(GFAsample(model,N))
  }
}

GFAsample <- function(model,N) {
  #
  # Generate data from a trained GFA model
  #
  # Inputs:
  #  model : a model trained by GFA()
  #  N     : the number of samples to be generated
  #
  # Outputs:
  # A list containing the following elements
  #  Y : A list of M elements, N times D[m] matrices
  #  Z : The latent variables used to generate the data, N times K
  #

  # Latent variables are white noise
  Z <- matrix(rnorm(N*model$K,0,1),N,model$K)
  
  Y <- vector("list",length=length(model$D))
  for(view in 1:length(model$D)) {
    # The mean is given by W times Z, and the noise is diagonal
    Y[[view]] <- Z %*% t(model$W[[view]]) +
      matrix(rnorm(N*model$D[view],0,1/sqrt(model$tau[view])),N,model$D[view])
  }

  return(list(Y=Y, Z=Z))
}
