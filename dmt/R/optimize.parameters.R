
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     

# Part of the inhumanity of the computer is that, once it is competently
# programmed and working smoothly, it is completely honest.
# - Isaac Asimov 


optimize.parameters <- function (X, Y, zDim = 1, priors = NULL, 
                                 marginalCovariances = "full", 
				 epsilon = 1e-6, convergence.steps = 3, verbose = FALSE) {

  # Suitable for at least:
  # nonmatched, prior$W, full marginals

  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) but there is option to set prior on W
  # (W.prior)

  # convergence.steps: convergence criteria need to be met at at least
  # 		       this many consecutive iteration steps.

  if ( verbose ) { cat("Initialize\n") }
  inits <- initialize2(X, Y, zDim, marginalCovariances)
  phi <- inits$phi
  phi.inv <- inits$phi.inv
  W <- inits$W
  Dcov <- inits$Dcov
  Dim <- inits$Dim
  nullmat <- inits$nullmat
  Nsamples <- inits$Nsamples

  # FIXME: handle priors completely outside this function later!

  if ( length(priors) == 0 ) { priors <- list() }

  ###  Wx ~ Wy prior inits  ###

  if ( verbose ) { cat("Checking the priors\n") }
  
  if ( !is.null(priors$Nm.wxwy.mean) ) {
    if ( length(priors$Nm.wxwy.mean) == 1 ) { priors$Nm.wxwy.mean <- priors$Nm.wxwy.mean * diag(1, nrow(X), nrow(Y)) }
    if ( ncol(priors$Nm.wxwy.mean) != nrow(X)){ stop("columns of priors$Nm.wxwy.mean must match rows of X") }
    if ( nrow(priors$Nm.wxwy.mean) != nrow(Y)){ stop("rows of priors$Nm.wxwy.mean must match rows of Y") }  
  }

  if ( is.null(priors$Nm.wxwy.sigma) || priors$Nm.wxwy.sigma == Inf ) {
    # Wx, Wy relation not constrained
    if ( verbose ) { cat("Wx ~ Wy free\n") }
    priors$Nm.wxwy.sigma <- Inf 
    # cost.W.exponential accepts also priors$W = NULL i.e. no W prior

    cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)

  } else if (priors$Nm.wxwy.sigma > 0) { # Wx ~ Wy constrained

    if ( verbose ) {cat("Wx ~ Wy constrained\n")}

    priors$T.tmp <- 1/(2 * Nsamples * priors$Nm.wxwy.sigma)
    # We assume here that Wy = T%*%Wx. Optimizing also T.
    T <- array(rnorm(Dim$Y*Dim$X,0,sd(W$X)), dim = c(Dim$Y, Dim$X))
    # Ensure that Wy = T * Wx:
    W$Y <- T%*%W$X
    # W not necessarily positive here

    cost.new <- cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov)

  } else if (priors$Nm.wxwy.sigma == 0) { # Wx = Wy

    if ( verbose ) { cat("Wx = Wy \n") }

    # Ensure that the dimensionality of given w matches with given zDim
    w <- as.matrix(inits$W$X[, 1:zDim], ncol = zDim)
    W <- list(X = w, Y = w, total = rbind(w, w))
    if ( !is.null(priors$W) ) {
     if ( verbose ) { cat(paste("prior for W: ", priors$W, "\n")) }
      cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)    
    } else {
      if ( verbose ) { cat(paste("no prior for W. \n")) }
      cost.new <- cost7(as.vector(W$X), phi, Dcov, Dim, priors)        
    }  
  }
  
  ###################################################

  # FIXME: remove par.change from function input

  par.changes <- rep(1e300, convergence.steps)

  if ( verbose ) { cat(paste("Starting iterations \n")) }

  # Convergence ends when consecutive cost function changes are below
  # the threshold and the last step reduces the cost
  while (any(par.changes > epsilon) || (par.changes[[convergence.steps]] < 0)) {

    if ( verbose ) { cat(cost.new); cat("\n") }

    cost.old <- cost.new

    ###################################################

    # Update W: initialize with previous W	

    W.old <- W

    if (!is.null(priors$W)) {
      # Regularized W

      if ( is.null(priors$Nm.wxwy.sigma) || priors$Nm.wxwy.sigma == Inf ) {

        # optimizes Wx and Wy assuming they are independent
        opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W.exponential, 
	             method = "L-BFGS-B", phi = phi, priors = priors, 
		     Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), 
		     lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

        # Convert optimized W parameter vector to actual matrices
        # Note that here we always assume that W is positive
        W <- get.W.nonneg(opt$par, Dim)
	
      } else if ( priors$Nm.wxwy.sigma == 0 ) {
      
      	# assuming Wx = Wy, we can speed up (FIXME; use analytical alternatives?)

        # SimCCA Wx = Wy with regularized W (W>=0)
        # message("Case Wx = Wy and regularized W.")

	opt <- optim(as.vector(W$X), cost7, method = "L-BFGS-B", 
	             phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, 
		     control = list(maxit = 1e6), 
                     lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

        W$X <- W$Y <- get.W.nonneg.identical(opt$par, Dim)
	W$total <- rbind(W$X, W$Y)
		
      } else {
        stop("W regularization implemented only for identical or independent Wx, Wy ie. priors$Nm.wxwy.sigma = 0 and priors$Nm.wxwy.sigma = Inf")
        # FIXME add the intermediates, should be straightforward by combining penalized optimizations
      }

    } else if (is.null(priors$W)) { # Unconstrained W

      if ( priors$Nm.wxwy.sigma == 0 ) { # Wx = Wy

        # assuming Wx = Wy
        # see Bach-Jordan 2005, sec. 4.1 for details
	# equations modified from there to match Wx = Wy case

        W <- W.simcca.EM(W, phi, Dim, Dcov)

      } else if ( priors$Nm.wxwy.sigma > 0 && priors$Nm.wxwy.sigma < Inf ) { # Wx ~ Wy constrained

        # Update W: initialize with previous W			       
        opt <- optim(c(as.vector(W$X), as.vector(T)), cost.W, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov,
               control = list(maxit = 1e6),lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
      
        # Convert optimized W parameter vector to actual matrices
        wt <- get.W(opt$par, Dim)
        W <- wt$W
        T <- wt$T		
      } else { # Wx, Wy, independent a priori -> pCCA
        stop("Special case, corresponding to pCCA. No need to loop over W, phi in optimization. Use pCCA function for direct solution.")
      }
    }
        
    W.new <- W # redundant?

    ##################################################

    # Update phi

    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)    

    if (marginalCovariances == "full") {

      phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

      if ( priors$Nm.wxwy.sigma == 0 ) { # Wx = Wy

        # FIXME: implement this also for other covariance structures

        # see Bach-Jordan 2005, sec. 4.1 for details
        M <- solve(t(W.old$X)%*%(phi.inv$X + phi.inv$Y)%*%W.old$X + diag(zDim))
        # beta <- M%*%t(W$X)%*%phi.inv.sum # ct. set.beta.fullcov(M, W$total, phi.inv$total)
        # FIXME: replace this with the general M.set functions
        # FIXME: speedup by sharing M/beta with W.simcca.EM?

        # Update phi
        phi <- phi.EM.simcca(Dcov, W.new, phi.inv, W.old, M)
	     
      } else {  # assuming in general Wx != Wy
      
        # also check from optimize.fullcov.R
        M <- set.M.full2(W.old, phi.inv) 
        phi <- phi.EM.cca(Dcov, W.new, phi.inv, W.old, M, nullmat)
      }

   } else if (marginalCovariances == "isotropic") {

        # M and beta Possibly useful for speedups later, when 
	# considering joint analysis of W and phi
        # set.M is for isotropic X or Y
	# these should OK without modifications here.
        #M <- list()
        #M$X <- set.M(W$X, phi$X)
        #M$Y <- set.M(W$Y, phi$Y)      
        #beta <- list()
        #beta$X <- set.beta(M$X, W$X, phi$X)
        #beta$Y <- set.beta(M$Y, W$Y, phi$Y)
        # FIXME: if this is not used, remove M and beta
        # however check their use in W update
        #phi <- update.phi(Dcov, M, beta, W, phi)
	
	# FIXME: speed up by using scalars as in the ppca/pfa/pcca
        phi.scalar.x <- unique(diag(phi$X))
        phi.scalar.y <- unique(diag(phi$Y))	
        phi$X <- update.phi.isotropic(Dcov$X, W$X, phi.scalar.x, Dim$X)
        phi$Y <- update.phi.isotropic(Dcov$Y, W$Y, phi.scalar.y, Dim$Y)

        # convert to matrices
	phi$X <- diag(phi$X, Dim$X)
	phi$Y <- diag(phi$Y, Dim$Y)
	phi$total <- diag(c(diag(phi$X), diag(phi$Y)))

        # FIXME: update.phi.isotropic possibly also usable here
	# but perhaps slower; update.phi is just an empirical estimate
	# the ML used in identical isotropic would be better

     } else if (marginalCovariances == "identical isotropic") {

        # FIXME: perhaps using similar implementation with "isotropic"
	# (separately for X and Y)
	# would be faster and about as accurate? test.

        # Calling set.M.isotropic
	# FIXME: speed up by using scalars as in the ppca/pfa/pcca
        phi.scalar <- unique(diag(phi$total))
        phi.estimate <- update.phi.isotropic(Dcov$total, W$total, phi.scalar, Dim$X + Dim$Y)
	#phi <- list(X = phi.estimate, Y = phi.estimate)

        # convert to matrices
	phi$X <- diag(phi.estimate, Dim$X)
	phi$Y <- diag(phi.estimate, Dim$Y)
	phi$total <- diag(phi.estimate, Dim$X + Dim$Y)

        # FIXME could be sped up by using scalars here, and similar treatment with W updates than with the "isotropic" option
	    
     } else if ( marginalCovariances == "diagonal" ) {

       phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

       # FIXME: speedups possible when Wx = Wy. Implement.     
       # FIXME: compare M, beta to isotropic/full cases and join common parts
       # FIXME needs to be checked!

       phi <- phi.diagonal.double(W$total, phi.inv$total, Dcov$total, Dim) 
       #phi$X <- phi.diagonal.single(W$total, phi.inv$total, Dcov$X, Dim)       
       #phi$Y <- phi.diagonal.single(W$total, phi.inv$total, Dcov$Y, Dim)            

     } else {
       stop("Unknown marginalCovariances parameter!")
     }

    ##########################################################################

    # MONITORING CONVERGENCE

    if (!is.null(priors$W)) { # W regularized

      # FIXME: cost7 and cost.W.exponential should both optimize nonneg W; combine?
      if (priors$Nm.wxwy.sigma == 0) {
        cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)
      } else if (priors$Nm.wxwy.sigma == Inf) {
        cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)      
      } else {
        stop("W regularization implemented only for independent and identical Wx, Wy cases.")
        # FIXME add the intermediates also
      }
    
    } else { # W not regularized; Wx ~ Wy is regularized

      if (priors$Nm.wxwy.sigma == 0) { # Extreme case Wx = Wy
        #message("Corresponds to simcca: W free; Wx = Wy.")
        cost.new <- cost7(as.vector(W$X), phi, Dcov, Dim, priors)    	
      } else if (priors$Nm.wxwy.sigma > 0 && priors$Nm.wxwy.sigma < Inf) {
        cost.new <- cost.W(c(as.vector(W$X), as.vector(T)), phi, priors, Dim, Dcov)
      } else if  (priors$Nm.wxwy.sigma == Inf)  {
      	stop("Corresponds to pCCA. No need for (slow) iterative optimization. Use pcca function instead.")
      } else {
        stop("Provide proper (nonneg. real) value for priors$Nm.wxwy.sigma")
      }
    }

    # remove the first element, add the new cost change into end
    # this way we keep track of the last congergence.steps iterations
    par.changes <- c(par.changes, (cost.old - cost.new))[-1]

  }

  if ( verbose ) {
    cat("par.changes")
    cat(par.changes)
    cat('\n')
    cat(paste("Iterations OK.\n"))
    
  }

  # FIXME
  # Needed later if phis are treated as scalars
  #if (marginalCovariances == "isotropic" || marginalCovariances == "identical isotropic") {
  #  # force these scalars into diagonal matrices
  #  phi$X <- diag(phi$X, nrow(X))
  #  phi$Y <- diag(phi$Y, nrow(Y))
  #  phi$total <- diag(c(diag(phi$X),diag(phi$Y)),(nrow(X)+nrow(Y)))
  #}

  W$total <- rbind(W$X, W$Y)  
  rownames(W$X) <- rownames(X)
  rownames(W$Y) <- rownames(Y)
  rownames(W$total) <- c(rownames(X), rownames(Y))
  
  rownames(phi$X) <- colnames(phi$X) <- rownames(X)
  rownames(phi$Y) <- colnames(phi$Y) <- rownames(Y)
  rownames(phi$total) <- colnames(phi$total) <- c(rownames(X), rownames(Y))

  return( list(W = W, phi = phi) )

}











