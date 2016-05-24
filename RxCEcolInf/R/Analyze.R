"form.mcmc.object" <-
  function(holder, names, title, burnin, mcmc, thin, ...) {
#    cat("burnin = ", burnin, "\n")
#    cat("mcmc = ", mcmc, "\n")
#    cat("thin = ", thin, "\n")
#    cat("nrow(holder) = ", nrow(holder), "\n")
    
    output <- mcmc(data=holder, start=(burnin+1),
                   end=mcmc,
                   thin=thin)
    varnames(output) <- as.list(names)
    attr(output,"title") <- title
    
    attribs <- list(...)
    K <- length(attribs)
    attrib.names <- names(attribs)

    if (K>0){
      for (i in 1:K){
        attr(output, attrib.names[i]) <- attribs[[i]]
      }
    }
    
    return(output)  
  }




`Analyze` <-
  function(fstring, rho.vec, data=NULL,
           num.iters=1000000, save.every=1000, burnin=10000,
           mu.vec.0=rep(log((.45/(mu.dim-1))/.55), mu.dim),
           kappa=10, nu=(mu.dim+6), psi=mu.dim,
           mu.vec.cu=runif(mu.dim, -3, 0), NNs.start=NULL,
           THETAS.start=NULL, prob.re=0.15, 
           sr.probs=NULL, sr.reps=NULL, keep.restart.info=FALSE,
           keepNNinternals=0, keepTHETAS=0,  
           nolocalmode=50, numscans=1, Diri=100, dof=4, 
           print.every=10000, debug=1)
{


  data.mats <- get.mat(fstring, data)
  ##print(colnames(data.mats$rowmat))
  ##print(colnames(data.mats$colmat))
  if (is.null(colnames(data.mats$rowmat))){
    rownames.local <- paste("r", 1:ncol(data.mats$rowmat), sep="")
  }
  else{
    rownames.local <- colnames(data.mats$rowmat)
  }
  if (is.null(colnames(data.mats$rowmat))){
    colnames.local <- paste("c", 1:ncol(data.mats$colmat), sep="")    
  }
  else{
    colnames.local <- colnames(data.mats$colmat)
  }  
  NNtots <- cbind(data.mats$rowmat, data.mats$colmat)
  nrow.pt <- ncol(data.mats$rowmat)
  ncol.pt <- ncol(data.mats$colmat)
  mu.dim <- (ncol.pt - 1) * nrow.pt
  ntables <- nrow(data.mats$colmat)
  
  # Check there are at least 2 rows/cols:
  if (nrow.pt<2 || ncol.pt<2)
    stop("Table must be at least 2x2")
  
  # Check the NNtots is valid:
  if (any(apply(NNtots[,1:nrow.pt],1,sum) != apply(NNtots[,(nrow.pt+1):(nrow.pt+ncol.pt)],1,sum)))
    stop("Row and column totals do not sum to the same total total")

  # Check valid iteration specs:
  if (burnin>num.iters) # NOTE: THIS IS NOT AN ADEQUATE ERROR CHECK
    stop("burnin must not be larger than num.iters")

#  cat("mu.dim = ", mu.dim, "\n")
#  print(mu.vec.0)
  
# Added 12/01/08:
# Temporary patch-up:
NNtots.list <- list(NULL)
mu.vec.cu.list <- list(NULL)
for (i in 1){
  NNtots.list[[i]] <- NNtots
  mu.vec.cu.list[[i]] <- mu.vec.cu
}

		# Extract the number of temperature rungs:
    N <- length(NNtots.list)

		#
    # Pre-computation:
		#

    # NNtots.list and mu.vec.cu.list should be lists for each temp rung:
		# NNbounds is createdas a list with bounds for each temperature rung:
    NNbounds.list <- list(NULL)
		# Error check mu.vec.cu list:
    if ((!is.list(mu.vec.cu.list))||(length(mu.vec.cu.list)!=N))
      stop("mu.vec.cu.list must be a list of the same length as NNtots.list")
    mu.len <- nrow.pt*(ncol.pt-1)
    if (length(mu.vec.0)!=mu.len)
      stop("mu.vec.0 must be of length nrow.pt")
		# Loop over temperature rungs:
    for (i in 1:N){
			# Error check mu.vec.cu list elements:
      if (length(mu.vec.cu.list[[i]])!=mu.len)
        stop("mu.vec.cu.list must have all elements of length nrow.pt")
			# Compute bounds for rung i:
      NNbounds.list[[i]] <- gq.bounds(nrow.pt,ncol.pt,NNtots.list[[i]])
			# Force correct storage types:
      storage.mode(NNtots.list[[i]]) <- "double"
      storage.mode(NNbounds.list[[i]]) <- "double"
      storage.mode(mu.vec.cu.list[[i]]) <- "double"
      NNtots.list[[i]] <- as.matrix(NNtots.list[[i]])
      NNbounds.list[[i]] <- as.matrix(NNbounds.list[[i]])
      mu.vec.cu.list[[i]] <- as.matrix(mu.vec.cu.list[[i]])
    }

		# Function to determine the choice of sr.reps:
		# x<0.11 => nchg
		# 0.11<=x <0.20 => 1 rep multinomial
		# 0.20<=x=<0.80 => 2 rep multinomial
		# x>0.80 => nchg
		#
		rep.func <- function(x){
			r <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
			for (i in 1:nrow(x)){
				for (j in 1:ncol(x)){
					tmp <- ifelse(x[i,j]<0.20,2,1)
					tmp <- ifelse((x[i,j]<0.11)||(x[i,j]>0.80),0,tmp)
					r[i,j] <- tmp
				}
			}
			return(r)
		}

		if (any(is.null(sr.probs))){
			# Default according to row fractions:
			rfracs <- NNtots.list[[1]][,1:nrow.pt]
			rfracs <- rfracs/matrix(rep(apply(rfracs,1,sum),nrow.pt),ncol=nrow.pt)
			sr.probs <- rfracs
		}
		storage.mode(sr.probs) <- "double"
		if (any(is.null(sr.reps))){
			# Default according to rep.func (above):
			sr.reps <- rep.func(sr.probs)
		}
		storage.mode(sr.reps) <- "integer"

    # Is there a pre-specified starting state of NNs or random start?
    if (!(is.null(NNs.start))){
      for (i in 1:N){
				# User-supplied starting states. Error check them for each rung:
        storage.mode(NNs.start[[i]]) <- "double"
        NNs.start[[i]] <- as.matrix(NNs.start[[i]])
        if (!(nrow(NNs.start[[i]])==nrow(NNtots.list[[i]])))
          stop("Invalid number of rows in NNs.start")
        if (!(ncol(NNs.start[[i]])==(nrow.pt*ncol.pt)))
          stop("Invalid number of cols in NNs.start")
      }
      # Pre-specified starting state for NNs:
      rstart.NNs <- FALSE
    } else {
      # Random starting state for NNs:
      rstart.NNs <- TRUE
    }
    
    # Is there a pre-specified starting state of the THETAS or random start?
    if (!(is.null(THETAS.start))){
      for (i in 1:N){
				# User-supplied starting values. Error check each rung:
        storage.mode(THETAS.start[[i]]) <- "double"
        THETAS.start[[i]] <- as.matrix(THETAS.start[[i]])
        if (!(nrow(THETAS.start[[i]])==nrow(NNtots[[i]])))
          stop("Invalid number of rows in THETAS.start")
        if (!(ncol(THETAS.start[[i]])==(nrow.pt*ncol.pt)))
          stop("Invalid number of cols in THETAS.start")
      }
      # Pre-specified starting state for NNs:
      rstart.THETAS <- FALSE
    } else {
      # Random starting state for NNs:
      rstart.THETAS <- TRUE
    }

    # Error-check Diri:
    Diri <- rep(as.numeric(Diri),nrow(NNtots.list[[1]]))
    if (length(Diri)!=nrow(NNtots.list[[1]]))
      stop("Diri should be a scalar quantity")
    if (any(Diri<0.0))
      stop("Diri should be non-negative")
    storage.mode(Diri) <- "double"

    # Other miscellaneuous error checks:
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
    if ((prob.re<0)||(prob.re>1))
      stop("prob.re must be in [0,1]")

    args <- list("NNtots" = NNtots.list,
                 "NNbounds" = NNbounds.list,
                 "numrows_pt" = as.integer(nrow.pt),
                 "numcols_pt" = as.integer(ncol.pt),
                 "sr_probs" = as.matrix(sr.probs),
                 "sr_reps" = as.matrix(sr.reps),
                 "num_iters" = as.integer(num.iters),
                 "psi_0" = as.numeric(psi),
                 "kappa_0" = as.numeric(kappa),
                 "nu_0" = as.numeric(nu),
                 "nolocalmode" = as.numeric(nolocalmode),
                 "dof" = as.numeric(dof),
                 "mu_vec_0" = as.numeric(mu.vec.0),
                 "mu_vec_cu" = mu.vec.cu.list,
                 "prob_re" = as.numeric(prob.re),
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
    gq.output <- .Call("Analyze", PACKAGE="RxCEcolInf", args)

    if (keep.restart.info){
      # Save the state of the RNG in case gq.restart is needed:
      gq.output$restart.info$random.seed <- .Random.seed
    }

  output.names <- NULL
  mu <- gq.output$draws$mu
  mu.names <- paste("mu", rep(rownames.local, each=(length(colnames.local)-1)),
                    rep(colnames.local[-length(colnames.local)],
                        length(rownames.local)), sep=".")
  mu.names <- paste(mu.names, 1:length(mu.names), sep=".")
  output.names <- c(output.names, mu.names)
  
  Sigma <- gq.output$draws$Sigma
  Sigma.names <- colnames(Sigma)
  Sigma.names <- gsub("_", ".", Sigma.names)
  holder <- paste("sd", rep(rownames.local, each=(length(colnames.local)-1)),
                  rep(colnames.local[-length(colnames.local)],
                      length(rownames.local)), sep=".")
  Sigma.names[1:length(holder)] <- paste(holder, 1:length(holder), sep=".")
  output.names <- c(output.names, Sigma.names)

  NNs <- gq.output$draws$NNs
  NN.names <- paste("NN", rep(rownames.local, each=length(colnames.local)),
                    rep(colnames.local, length(rownames.local)), sep=".")
  output.names <- c(output.names, NN.names)
 # print(colnames(NNs))
 # print(NN.names)
  
  LAMBDA <- gq.output$draws$LAMBDA
  LAMBDA.names <- paste("LAMBDA",
                        rep(rownames.local,
                            each=(length(colnames.local) -1)),
                        rep(colnames.local[1:(length(colnames.local) - 1)],
                            length(rownames.local)),
                        sep=".")
  output.names <- c(output.names, LAMBDA.names)
#  print(colnames(LAMBDA))
#  print(LAMBDA.names)
  
  TURNOUT <- gq.output$draws$TURNOUT
  TURNOUT.names <- paste("TURNOUT", rownames.local, sep=".")
  output.names <- c(output.names, TURNOUT.names)
#  print(colnames(TURNOUT))
#  print(TURNOUT.names)
  
  
  GAMMA <- gq.output$draws$GAMMA
  GAMMA.names <- paste("GAMMA", rownames.local, sep=".")
  output.names <- c(output.names, GAMMA.names)
#  print(colnames(GAMMA))
#  print(GAMMA.names)

  BETA <- gq.output$draws$BETA
  BETA.names <- paste("BETA",
                        rep(rownames.local,
                            each=(length(colnames.local) - 1)),
                        rep(colnames.local[1:(length(colnames.local) - 1)],
                            length(rownames.local)),
                        sep=".")
  output.names <- c(output.names, BETA.names)
#  print(colnames(BETA))
#  print(BETA.names)

  outmat <- cbind(mu, Sigma, NNs, LAMBDA, TURNOUT, GAMMA, BETA)
  num.to.drop <- as.integer(burnin / save.every)
  if (num.to.drop !=0) outmat <- outmat[-c(1:num.to.drop),]
  
##  print(dim(outmat))
##  print(length(output.names))
  if (!is.null(gq.output$THETAS.saved)) {
    THETA <- t(gq.output$THETAS.saved)
    THETA.names <- paste(
                         rep(rownames.local,
                             each=(length(colnames.local))),
                         rep(colnames.local,
                             length(rownames.local)),
                         sep=".")
    THETA.names <- paste("THETA.table",
                         rep(1:ntables, each=(nrow.pt * ncol.pt)),
                         rep(THETA.names, ntables))
    THETA.thin <- as.integer(num.iters / nrow(THETA))
    num.to.drop <- as.integer(burnin / THETA.thin)
    if (num.to.drop != 0) THETA <- THETA[-c(1:num.to.drop),]
    THETA <- mcmc(data=THETA, start=(burnin+1), end=num.iters,
                  thin=THETA.thin)
    varnames(THETA) <- THETA.names
  } else {
    THETA <- NULL
  }
 
  if (!is.null(gq.output$NN.internals.saved)) {
    NN <- t(gq.output$NN.internals.saved)
    NN.names <- paste(
                      rep(rownames.local,
                            each=(length(colnames.local))),
                      rep(colnames.local,
                            length(rownames.local)),
                      sep=".")
    NN.names <- paste("NN.table",
                    rep(1:ntables, each=(nrow.pt * ncol.pt)),
                    rep(NN.names, ntables))
    NN.thin <- as.integer(num.iters / nrow(NN))
    num.to.drop <- as.integer(burnin / NN.thin)
    if (num.to.drop != 0) NN <- NN[-c(1:num.to.drop),]
    NN <- mcmc(data=NN, start=(burnin+1), end=num.iters,
             thin=NN.thin)
    varnames(NN) <- NN.names
  } else {
    NN <- NULL
  }
 


  
  mcmc.out <- form.mcmc.object(outmat, output.names, "RxCEcolInf Output",
                               burnin=burnin, mcmc=num.iters, thin=save.every,
                               acc.t=gq.output$acc.t,
                               acc.Diri=gq.output$acc.Diri,
                               vld.multinom=gq.output$vld.multinom,
                               acc.multinom=gq.output$acc.multinom,
                               numrows.pt=gq.output$numrows.pt,
                               numcols.pt=gq.output$numcols.pt,
                               THETA=THETA,
                               NN.internals=NN,
                               restart.info=gq.output$restart.info)
                               

  
  
  return(mcmc.out)
}
