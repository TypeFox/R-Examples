# (C) 2008-2011 Olli-Pekka Huovilainen and Leo Lahti 
# All rights reserved.
# FreeBSD license (keep this notice).



# "Computer Science is no more about computers than astronomy is about
#  telescopes." 
# - E. W. Dijkstra



fit.dependency.model <- function (X, Y,
          zDimension = 1,
          marginalCovariances = "full",
          epsilon = 1e-3,
          priors = list(), matched = TRUE,
          includeData = TRUE, calculateZ = TRUE, verbose = FALSE)
{

  # zDimension = 1; marginalCovariances = "full"; H = 1; sigmas = 0; epsilon = 1e-3; mySeed = 123; priors = NULL
  
  # Fits the generative model
  # X = Wx * z + epsx
  # Y = Wy * z + epsy
  # with various modeling assumptions

  if ( verbose ) { cat("Checking data\n") }
  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
  intercept <- dat$intercept
       
  # FIXME store/return intercepts as well; further dependency models including intercepts
  if ( nrow(X) < nrow(Y) ) {stop("If the two data matrices do not have equal dimensionality, place the smaller one in Y.")} # FIXME automate
  if ( !nrow(X) == nrow(Y) ) {
    #message("The data sets have unequal dimensionality.")
    if ( matched ) { stop("Cannot use matched methods for nonmatched data.") }
  }

  if ( verbose ) { cat("Checking inputs\n") }
  if ( !is.null(priors$Nm.wxwy.sigma ) && priors$Nm.wxwy.sigma == Inf) { matched <- FALSE; message("priors$Nm.wxwy.sigma == Inf; Wx ~ Wy independendent i.e. matched = FALSE") }  
  if ( epsilon == 0 )  { epsilon <- 1e-3 } # avoid numerical overflows
  res <- NA; method <- ""
  
  if (!matched) {

    if (verbose) { cat("Model for non-matched case\n") }

    if ( is.null(priors$Nm.wxwy.sigma )) {
      #warning("priors$Nm.wxwy.sigma not implemented for non-matched variables. Setting priors$Nm.wxwy.sigma = Inf.")
      priors$Nm.wxwy.sigma <- Inf
    }

    if ( is.null(priors$W) ) { 

      if ( verbose ) { cat("Wx ~ Wy free. No regularization for W.\n") }
      if ( verbose ) { cat(marginalCovariances); cat("\n") }
      
      if ( marginalCovariances == "full" ) { # standard pCCA

        method <- "pCCA"
        res <- calc.pcca(X, Y, zDimension)

      } else if (marginalCovariances == "diagonal") { 

        # pFA for two data sets corresponds to 
        # standard pFA for concatenated data

        res <- calc.pfa(X, Y, zDimension)    
        method <- "pFA"

      } else if (marginalCovariances == "isotropic") {
      
        # phiX != phiY in general
        # FIXME: add tests

        res <- calc.pcca.with.isotropic.margins(X, Y, zDimension)

        # force these scalars into diagonal matrices                  
        res$phi$X <- diag(res$phi$X, nrow(X))
        res$phi$Y <- diag(res$phi$Y, nrow(Y))           
        res$phi$total <- diag(c(diag(res$phi$X),diag(res$phi$Y)), (nrow(X) + nrow(Y)))

        # Update Phi                                         
	#phi.scalar.x <- unique(diag(phi$X))
	#phi.scalar.y <- unique(diag(phi$Y))	
        # modified from calc.pcca.with.isotropic.margins
        #phi$X <- update.phi.isotropic(Dcov$X, W$X, phi.scalar.x, Dim$X)
        #phi$Y <- update.phi.isotropic(Dcov$Y, W$Y, phi.scalar.y, Dim$Y) 
                                                     							      						       							
        method <- "pCCA with isotropic margins"
	
      } else if (marginalCovariances == "identical isotropic") {
      
        # FIXME: add tests for this
        # pPCA for two data sets corresponds to 
        # standard pPCA for concatenated data

        res <- calc.ppca(X, Y, zDimension)
        method <- "pPCA"
	
      } else { stop("Erroneous marginalCovariances parameter provided!") }

    } else if ( !is.null(priors$W) ) { 
      
      if ( verbose ) { cat("Wx ~ Wy free; exponential (nonnegative) prior for W.\n") }

      # Prior for W is given -> no analytical solution to EM
      # Exponential prior for W,
      # Used to force positive W with exponential distribution.
      # priors$W is the rate parameter of the exponential. 
      # The smaller, the flatter tail.

      # TODO: implement sparsity prior W ~ N(0, sd*I)

      res <- optimize.parameters(X, Y, zDim = zDimension, priors = priors,                                                                 
                                   marginalCovariances = marginalCovariances,           
                                   epsilon = epsilon, convergence.steps = 3, verbose = verbose)

      method <- "Free Wx ~ Wy with exponential priors for Wx and Wy. Check marginal covariances from parameters."
	
    }
    
  } else if (matched) {

    if ( verbose ) { cat("Matched features case\n") }
      
    # Matrix normal distribution variance not specified
    if ( is.null(priors$Nm.wxwy.sigma) ) {
      message("Matched variables but priors$Nm.wxwy.sigma not given, using strong matching with Wx = Wy.")
      priors$Nm.wxwy.sigma <- 0
    }
      
    # Matrix normal distribution mean matrix not specified
    if ( is.null(priors$Nm.wxwy.mean) || is.na(priors$Nm.wxwy.mean)) {
      message("The matrix Nm.wxwy.mean is not specified. Using identity matrix.")
      priors$Nm.wxwy.mean <- 1
    }    

    method <- "pSimCCA"
        
    # Case IIa: fully constrained case Wx = Wy
    if ( priors$Nm.wxwy.sigma == 0 ) { #Wx = Wy        
        
      if ( verbose ) { cat("Assuming Wx = Wy\n") }
	
      #  SimCCA with full covariances with constraint Wx = Wy
      #  "probsimcca.full" = aucs.simcca.full      
      #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
      #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
      #  number of features in X and Y)

      # If prior for W is given, we must optimize W (no analytical
      # solution to EM)
          
      # Regularization for W (W > 0 implemented)
      if (!is.null(priors$W)) {
            
	if ( verbose ) { cat("Wx = Wy with regularized W (W>=0)\n") }
	if ( verbose ) { cat(marginalCovariances); cat("\n") }	

        res <- optimize.parameters(X, Y, zDim = zDimension, priors = priors, 
	       			   marginalCovariances = marginalCovariances, 
				   epsilon = epsilon, convergence.steps = 3, verbose = verbose)

        method <- "pCCA with W prior"

      } else if (is.null(priors$W)) {
        
	if ( verbose ) { cat("Wx = Wy; free W.\n") }

          # mlsp'09 simcca
          # message("Case Wx = Wy. No regularization for W.")
	  
	  # use this for full W (EM algorithm, unstable for n ~ p)
         res <- optimize.parameters(X, Y, zDim = zDimension, priors = priors, 
                                   marginalCovariances = marginalCovariances,           
                                   epsilon = epsilon, convergence.steps = 3,
                                   verbose = verbose)


          method <- "matched case Wx = Wy with unconstrained W. Check covariances from parameters."
          # FIXME: speeups possible here when Wx = Wy but not yet implemented with other than full covs
      }
      
    } else if ( priors$Nm.wxwy.sigma > 0 ) {
      # Case IIb: partially constrained Wx ~ Wy
                
      if ( verbose ) { cat("partially constrained Wx ~ Wy.\n") }
		
      if ( !is.null(priors$W) ) {
        if ( verbose ) {cat("regularized W.\n")}
        # FIXME: consider adding fast option with simply nonnegative W but no distributional priors
        stop("Not implemented regularization for W with partially constrained Wx ~ Wy.")
      } else if (is.null(priors$W)) {
        if ( verbose ) { cat("Partially constrained Wx ~ Wy. No regularization for W.\n") }
        if ( verbose ) { cat(marginalCovariances); cat("\n") }            		  
			  
        if ( marginalCovariances == 'isotropic' ) {
	  # message("SimCCA with isotropic covariances and regularized H (through sigmas).")
	
          # FIXME: consider later adding other covariance structures if needed?
	  # note that the computation is slow then          		

         res <- optimize.parameters(X, Y, zDim = zDimension, priors = priors,
	     			       marginalCovariances = marginalCovariances,                                   
				       epsilon = epsilon, convergence.steps = 3, verbose = verbose)
                                                                    
         method <- "constrained Wx~Wy with matrix normal distribution prior"


        } else if ( !marginalCovariances == 'isotropic' ) {
          stop("Only isotropic marginal covariances implemented with constrained Wx ~ Wy in the general case.")
        }
      } 
    }
  }


  if ( verbose ) { cat("Checking the model..\n") }

  # Test whether model exists for given arguments
  if ( any(is.na(unlist(res))) ) {
    stop("Error with model parameters.")
  } else {
    params <- list(marginalCovariances = marginalCovariances, Nm.wxwy.mean = priors$Nm.wxwy.mean, Nm.wxwy.sigma = priors$Nm.wxwy.sigma, zDimension = zDimension, epsilon = epsilon)
    score <- dependency.score( res )
  }
  
  if ( verbose ) {cat("Creating DependencyModel object..\n")}
  
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)	
  if ( includeData ) model@data <- list(X = X, Y = Y)
  if ( calculateZ ) model@z <- z.expectation(model, X, Y)
  
  if ( verbose ) { cat("fit.dependency.model OK.\n") }
  
  model
}

# FIME: consider whether we should keep this or remove			    
#pcca.isotropic <- function(X, Y, zDimension = NULL, matched = FALSE, epsilon = 1e-6, includeData = TRUE, calculateZ = TRUE){
#  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }#
#
#  fit.dependency.model(X,Y,zDimension,marginalCovariances = "isotropic", epsilon = 1e-6,
#                       includeData = includeData, calculateZ = calculateZ)          
#}                                                                      

pcca <- function (X, Y, zDimension = NULL, includeData = TRUE, calculateZ = TRUE) {

  # (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti
  # License: FreeBSD (keep this notice)

  # replaces: solve.CCA.full
  
  # If zDimension given, then
  # only pick zDimension first principal components
  # and estimate marginals accordingly
  # relies on the fact that the principal components
  # can be estimated consecutively in pCCA
  
  # Add here centering of the data matrices X, Y
  # (center the dimensions to 0)

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension

  res <- calc.pcca(X, Y, zDimension)

  method <- "pCCA"
  params <- list(marginalCovariances = "full", zDimension = zDimension)
  score <- dependency.score( res )
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)
  if ( includeData ) model@data <- list(X = X, Y = Y)
  if ( calculateZ )  model@z <- z.expectation(model, X, Y) 
  model
  
}


calc.pcca <- function (X, Y, zDimension) {

  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))

  # Solve W (solve.w utilizes Archambeau06 equations from PCA for shortcut)
  # FIXME: compare accuracy and speed to direct EM update scheme?
  W <- solve.w(t(X), t(Y), Dcov$X, Dcov$Y, zDimension)

  # Then take only the zDimension first components if defined
  W$X <- as.matrix(W$X)
  W$Y <- as.matrix(W$Y)
  W$total <- rbind(W$X, W$Y)
        
  # estimate
  phi <- list()
  phi$X <- Dcov$X - W$X%*%t(W$X)
  phi$Y <- Dcov$Y - W$Y%*%t(W$Y)

  # Retrieve principal canonical components from the prob.CCA model
  # assuming that W and phi are known (at least in the full-rank case)
  # U <- solve.archambeau(X, Y, W$X, W$Y, phi$X, phi$Y)

  list(W = W, phi = phi)
  
}

ppca <- function (X, Y = NULL, zDimension = NULL, includeData = TRUE, calculateZ = TRUE) {

  dat <- check.data(X, Y, zDimension)
  X <- dat$X  
  Y <- dat$Y
  zDimension <- dat$zDimension

  res <- calc.ppca(X, Y, zDimension)
  
  method <- "pPCA"	
  params <- list(marginalCovariances = "isotropic", zDimension = zDimension)
  score <- dependency.score( res )
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)
  if (includeData) { model@data <- list(X = X, Y = Y) }
  if (calculateZ) { model@z <- z.expectation(model, X, Y) }
  model
 	    
}

# FIXME: now calc.ppca used only in one-data case by function ppca
# either remove two-data case, or test and compare with fit.dependency.model and take into use
calc.ppca <- function (X, Y, zDimension) {

  # Replaces function solve.CCA

  # if zDimension = NULL then full-rank solution is computed
  
  # Probabilistic PCA
  # (See Tipping and Bishop 1999)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # If only X is given in the argument, compute
  # pPCA for X

  # If both X and Y are given in the argument, compute
  # pPCA for concatenated [X; Y]
  # Assuming isotropic and identical marginal noise, the
  # principal subspace will capture the dependencies between X and Y.
  # This corresponds to the model (sigmax = sigmay = sigma)
  # X ~ N(Wx*z, sigma*I); Y ~ N(Wy*z, sigma*I)
  # This provides a simple comparison method for more
  # detailed dependency models.

  if ( is.null(Y) ) {
    res <- ppca.calculate(X, zDimension)
    phi <- list(total = diag(res$phi, nrow(X)))
    rownames(res$W) <- rownames(X)
    colnames(phi$total) <- rownames(phi$total) <- rownames(X)
    #W <- list(X = res$W, total = res$W)
    W <- list(total = res$W)    
  } else {
    # If second argument (Y) given, compute
    # pPCA with two concatenated data sets
    res <- ppca.calculate(rbind(X,Y), zDimension)

    # Make phi diagonal matrix
    phitotal <- diag(res$phi,(nrow(X) + nrow(Y)))

    # Variable names to W and phi
    rownames(res$W) <- c(rownames(X),rownames(Y))
    rownames(phitotal) <- c(rownames(X),rownames(Y))
    colnames(phitotal) <- rownames(phitotal)

    # Divide W and phi to X and Y parts
    phi <- list(X = phitotal[(1:nrow(X)),(1:nrow(X))], 
                Y = phitotal[-(1:nrow(X)),-(1:nrow(X))],
    	        total = phitotal)
    W <- list(X = as.matrix(res$W[(1:nrow(X)),]), 
              Y = as.matrix(res$W[-(1:nrow(X)),]), 
	      total = res$W)

  }
  # Note that if X, Y given then phi$X = phi$Y in the pCCA model
  # Here W corresponds to W$total of other functions when X, Y both given
  list(W = W, phi = phi)
}


ppca.calculate <- function (X, zDimension) {

 # FIXME: ensure that X is zero-mean

  # Probabilistic PCA
  # (See Tipping and Bishop 1999 / 3.2)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # Use full-rank if dimensionality is not specified
  zDimension <- ifelse(is.null(zDimension), nrow(X), zDimension)

  # eigenvalues D and eigenvectors U       
  duv <- svd(X)
  U <- duv$u                                                 
  D <- sort(duv$d, decreasing = TRUE)  

  # ML estimate for phi (in pPCA)
  phi <- sum(D[-seq(zDimension)])/(nrow(duv$u) - zDimension)

  # ML estimate for W, given phi (in pPCA)
  # Here set R <- I (R is an arbitrary orthogonal rotation matrix)  
  W <- as.matrix(U[, 1:zDimension])%*%sqrt(diag(D)[seq(zDimension), seq(zDimension)] - 
       		      phi * diag(1, zDimension, zDimension))

  if (zDimension == nrow(X)) {

    # If W is full-rank then the isotropic error term will disappear, assuming data X is gaussian
    # then X%*%t(X) i.e. cov(t(X)) approximates W%*%t(W) (since X ~ Wz and z ~ N(0,I))

    # Note rotational ambiguity for W, Z 
    cat("Full-rank PCA calculated, isotropic error term is zero.\n")
    W <- matrix.sqrt(cov(t(X)))
    phi <- 0
  }
  
  list(W = W, phi = phi)
}

pfa <- function (X, Y = NULL,
                 zDimension = NULL,
                 includeData = TRUE,
                 calculateZ = TRUE, priors = NULL) {

  # Probabilistic factorial analysis model as proposed in
  # EM Algorithms for ML Factoral Analysis, Rubin D. and 
  # Thayer D. 1982

  # Assumption: R = I

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
  method <- "pFA"
  params <- list(marginalCovariances = "diagonal", zDimension = zDimension)

  if (nrow(X) == 1 && (nrow(Y) == 1 || is.null(Y)) ) {

    score <- Inf
    res <- list(W = dat, phi = list(X = 0, Y = 0))
    if (is.null(Y)) {res$phi$Y <- res$phi$X <- NULL; res$phi$total <- 0}    

  } else {

    res <- calc.pfa(X, Y, zDimension, priors)
    score <- dependency.score( res )  

  }

  model <- new("DependencyModel",
               W = res$W, phi = res$phi,
               score = score,
               method = method,
               params = params)
  if ( includeData ) model@data <- list(X = X, Y = Y)

  if ( calculateZ ) {
    if (nrow(X) == 1 && (nrow(Y) == 1 || is.null(Y)) ) {
      model@z <- X
    } else {
      model@z <- z.expectation(model, X, Y)
    }
  }

  model

}

calc.pfa <- function (X, Y, zDimension, priors = NULL) {

  # Y.rubin is Y in (Rubin & Thayer, 1982)
  # Variables on columns and samples on rows
  # W corresponds to t(beta)
  if (is.null(Y)){
    Y.rubin <- t(X)
    # Factor loading matrix
    Wt <- t(eigen(cov(t(X)))$vectors[, 1:zDimension])
  } else {
    Y.rubin <- cbind(t(X), t(Y))
    # Use different initialization for Wt when data has inequal dimensionalities
    if (nrow(X) != nrow(Y)) {
      Wt <- t(eigen(cov(Y.rubin))$vectors[,1:zDimension])
    } else {
      init <- initialize2(X, Y, zDimension, marginalCovariances = "diagonal")
      # Factor loading matrix
      Wt <- t(init$W$total[,1:zDimension])
    }
  }
  
  epsilon <- 1e-3
  colnames(Wt) <- colnames(Y.rubin)
  tau2 <- diag(ncol(Y.rubin))

  Cyy <- cov(Y.rubin)
  delta <- 1e12
  # EM
  while(delta > epsilon){

    Wt.old <- Wt
    tau2.old <- tau2

    # E-step
    invtau2 <- solve(tau2)
    binv <- Wt%*%invtau2
    tbb <- invtau2 - (invtau2%*%t(Wt))%*%solve(diag(zDimension) + binv%*%t(Wt))%*%(binv)
    d <- tbb%*%t(Wt)
    D <- diag(zDimension) - Wt%*%d
    cyd   <- Cyy%*%d
    
    # M-step
    # Update W
    # FIXME: combine calc.pca, calc.pfa and calc.cca into one uniform model?
    if (is.null(priors)) {
      #message("Analytical optimization")      
      Wt  <- solve(t(d)%*%cyd + D)%*%t(cyd) # WORKS    
      # Also obtained with numerical optimization:
      # Wt <- update.W.singledata(Wt, X, tau2)
    } else if (!is.null(priors$W) && priors$W > 0) {
      #message("Numerical optimization")
      Wt <- update.W.singledata(Wt, X, tau2, priors)
    }
    
    # Update margin/s 
    tau2  <- diag(diag(Cyy - cyd%*%solve(t(d)%*%cyd + D)%*%t(cyd)))
    # Check cost function convergence
    delta <- max(sum(abs(tau2 - tau2.old)), sum(abs(Wt - Wt.old)))

  }

  # Convert names as same in other methods
  if ( is.null(Y) ){
      W <- list(total = t(Wt))
    phi <- list(total = tau2)
  } else {
      W <- list(X = as.matrix(t(Wt)[(1:nrow(X)),]), Y = as.matrix(t(Wt)[-(1:nrow(X)),]), total = t(Wt))
    phi <- list(X = tau2[1:nrow(X),1:nrow(X)], Y = tau2[-(1:nrow(X)),-(1:nrow(X))], total = tau2)                
  }
  
  list(W = W, phi = phi)

}




Pcca.with.isotropic.margins <- function (X, Y, zDimension = 1, epsilon = 1e-6, delta = 1e6) {

  # epsilon and delta are convergence parameters
  # zDimension determines the dimensionality of the shared latent variable Z

  #  Dependency model
  #  X ~ N(Wx*z, sigmax*I)
  #  y ~ N(Wy*z, sigmay*I)
  #  i.e. isotropic marginals but in general  sigmax != sigmay
  # This is a special case of probabilistic factor analysis model

  # FIXME: ensure that X, Y have zero-mean (shift if necessary);
  # alternatively add mean parameter in the model

  res <- calc.pcca.with.isotropic.margins(X, Y, zDimension, epsilon = epsilon, delta = delta)
  phi <- res$phi  
    W <- res$W   


  colnames(phi$X) <- rownames(phi$X) <- rownames(X)
  colnames(phi$Y) <- rownames(phi$Y) <- rownames(Y)
  colnames(phi$total) <- rownames(phi$total) <- c(rownames(X), rownames(Y))
  
  list(phi = phi, W = W)

  # FIXME provide here proper DependencyModel object as in pcca, pfa and ppca
}


calc.pcca.with.isotropic.margins <- function (X, Y, zDimension, epsilon = 1e-6, delta = 1e6) {

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
  
  # initialize
     inits <- initialize2(X, Y, zDimension, marginalCovariances = "isotropic")
      Dcov <- inits$Dcov
       Dim <- inits$Dim
         W <- inits$W  

  # FIXME: ensure that X, Y have zero-mean (shift if necessary);
  # alternatively add mean parameter in the model
  phi <- list(X = 1, Y = 1)
  
  # iterate until convergence:
  while (delta > epsilon) {

    W.old <- W
          
    ##########################################

    # Update Phi
    phi$X <- update.phi.isotropic(Dcov$X, W$X, phi$X, Dim$X) 
    phi$Y <- update.phi.isotropic(Dcov$Y, W$Y, phi$Y, Dim$Y)
          
    #######################################

    # Full CCA update for W 

        phi.inv <- list()
        phi.inv$X <- diag(rep(1/phi$X, Dim$X))
        phi.inv$Y <- diag(rep(1/phi$Y, Dim$Y))
        phi.inv$total <- diag(c(rep(1/phi$X, Dim$X), rep(1/phi$Y, Dim$Y)))

         #M <- set.M.full(W$total, phi.inv.full) # corresponds to G in Bishop's book
          M <- set.M.full2(W, phi.inv) # modified for G in Bishop's book	  
       beta <- set.beta.fullcov(M, W$total, phi.inv$total)
    W$total <- W.cca.EM(Dcov, M, beta)
        W$X <- matrix(W$total[1:Dim$X,], nrow = Dim$X)
        W$Y <- matrix(W$total[-(1:Dim$X),], nrow = Dim$Y)

    ########################################
          
    # check convergence (enough to check W)
    delta <- max(abs(as.vector(W$total - W.old$total)))
          
  }

  # FIXME: add 'intercept' field to DependencyModel?
  list(W = W, phi = phi)

}

