`gq.restart` <-
  function(gq.output,extra.iter,print.every=1000,keep.restart.info=TRUE,debug=1)
{
    if (extra.iter==0){
      return(gq.output)
    }
    # Pre-processing:
    gq.restart.list <- attr(gq.output,"restart.info")
    if (is.null(gq.restart.list))
      stop("No restart information provided, did you select keep.restart.info=TRUE in your original run?")

# Added 12/01/08:
# Temporary patch-up:
mu.vec.cu.list <- list(NULL)
for (i in 1){
  mu.vec.cu.list[[i]] <- gq.restart.list$mu.vec.cu
}
    # Extract information for the new run:
    NNtots.list   <- gq.restart.list$NNtots
    NNbounds.list <- gq.restart.list$NNbounds

    # Patch up:
    N <- 1 
		# Loop over temperature rungs:
    for (i in 1:N){
			# Force correct storage types:
      storage.mode(mu.vec.cu.list[[i]]) <- "double"
      storage.mode(NNtots.list[[i]]) <- "double"
      storage.mode(NNbounds.list[[i]]) <- "double"
      NNtots.list[[i]] <- as.matrix(NNtots.list[[i]])
    }
    nrow.pt <- attr(gq.output,"numrows.pt")
    ncol.pt <- attr(gq.output,"numcols.pt")
    num.iters <- extra.iter
    psi <- gq.restart.list$psi.0
    kappa <- gq.restart.list$kappa.0
    nu <- gq.restart.list$nu.0
    nolocalmode <- gq.restart.list$nolocalmode
    dof <- gq.restart.list$dof
    mu.vec.0 <- gq.restart.list$mu.vec.0
    Diri <- gq.restart.list$use.Diri.every.vec
    rho.vec <- gq.restart.list$rho.vec
    save.every <- gq.restart.list$save.every
    numscans <- gq.restart.list$numscans
    keepNNinternals <- gq.restart.list$keepNNinternals
    keepTHETAS <- gq.restart.list$keepTHETAS
    NNs.start <- gq.restart.list$NNs
    storage.mode(NNs.start) <- "double"
    # Not a random starting state for NNs:
    rstart.NNs <- FALSE
    THETAS.start <- gq.restart.list$THETAS
    storage.mode(THETAS.start) <- "double"
    # Not a random starting state for NNs:
    rstart.THETAS <- FALSE

    # Set the random seed:
    set.seed(gq.restart.list$random.seed)

    # No longer needed:
    gq.restart.list <- NULL

    # Error checks:
    if ((nrow.pt<0)||(ncol.pt<0))
      stop("Number of rows and columns must be positive integers")
    if (num.iters<=0)
      stop("Please specify a positive number of iterations")
    if (debug<0)
      dbg <- 0
    if (print.every<=0)
      stop("print.every must be positive integer")
    if (numscans<1)
      stop("numscans must be positive integer")

    args <- list("NNtots" = NNtots.list,
                 "NNbounds" = NNbounds.list,
                 "numrows_pt" = as.integer(nrow.pt),
                 "numcols_pt" = as.integer(ncol.pt),
                 "num_iters" = as.integer(num.iters),
                 "psi_0" = as.numeric(psi),
                 "kappa_0" = as.numeric(kappa),
                 "nu_0" = as.numeric(nu),
                 "nolocalmode" = as.numeric(nolocalmode),
                 "dof" = as.numeric(dof),
                 "mu_vec_0" = as.numeric(mu.vec.0),
                 "mu_vec_cu" = mu.vec.cu.list,
                 "use_Diri_every_vec" = as.numeric(Diri),
                 "rho_vec" = as.numeric(rho.vec),
                 "save_every" = as.numeric(save.every),
                 "numscans" = as.integer(numscans),
                 "rstart_NNs" = as.integer(rstart.NNs),
                 "NNs_start" = NNs.start,
                 "rstart_THETAS" = as.integer(rstart.THETAS),
                 "THETAS_start" = THETAS.start,
                 "keep_restart_info" = as.integer(keep.restart.info),
                 "keepNNinternals" = as.double(keepNNinternals),
                 "keepTHETAS" = as.double(keepTHETAS),
                 "print_every" = as.integer(print.every),
                 "dbg" = as.integer(debug))

    gq.append.output <- .Call("Analyze", PACKAGE="RxCEcolInf", args)

    # Mesh the two outputs together:
    attr(gq.output,"restart.info")  <- attr(gq.append.output,"restart.info")
    gq.output$draws$mu      <- rbind(gq.output$draws$mu,gq.append.output$draws$mu)
    gq.output$draws$Sigma   <- rbind(gq.output$draws$Sigma,gq.append.output$draws$Sigma)
    gq.output$draws$NNs     <- rbind(gq.output$draws$NNs,gq.append.output$draws$NNs)
    gq.output$draws$LAMBDA  <- rbind(gq.output$draws$LAMBDA,gq.append.output$draws$LAMBDA)
    gq.output$draws$TURNOUT <- rbind(gq.output$draws$TURNOUT,gq.append.output$draws$TURNOUT)
    gq.output$draws$GAMMA   <- rbind(gq.output$draws$GAMMA,gq.append.output$draws$GAMMA)
    if (keep.restart.info){
      # Save the state of the RNG in case gq.restart is needed:
      attr(gq.output,"restart.info")$random.seed <- .Random.seed
    }

    return(gq.output)
}
