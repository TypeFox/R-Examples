
`AnalyzeWithExitPoll` <-
  function(fstring,  rho.vec,  exitpoll, data=NULL,
           num.iters=1000000, save.every=1000, burnin=10000,
           mu.vec.0=rep(log((0.45/(mu.dim-1))/0.55), mu.dim),
           kappa=10, nu=(mu.dim+6), psi=mu.dim,
           mu.vec.cu=runif(mu.dim, -3, 0),  NNs.start=NULL,
           MMs.start=NULL, THETAS.start=NULL,
           sr.probs=NULL, sr.reps=NULL, keep.restart.info=FALSE,
           keepNNinternals=0, keepTHETAS=0, 
           nolocalmode=50, numscans=1, Diri=100, dof=4, eschew=FALSE,
           print.every=10000, debug=1)
{

  data.mats <- get.mat(fstring, data)
  ##print(colnames(data.mats$rowmat))
  ##print(colnames(data.mats$colmat))
  if (is.null(colnames(data.mats$rowmat))){
    rownames.local <- paste("r", 1:ncol(data.mats$rowmat), sep="")
  } else {
    rownames.local <- colnames(data.mats$rowmat)
  }
  if (is.null(colnames(data.mats$rowmat))){
    colnames.local <- paste("c", 1:ncol(data.mats$colmat), sep="")    
  } else {
    colnames.local <- colnames(data.mats$colmat)
  }  
  NNtots <- cbind(data.mats$rowmat, data.mats$colmat)
  nrow.pt <- ncol(data.mats$rowmat)
  ncol.pt <- ncol(data.mats$colmat)
  mu.dim <- (ncol.pt - 1) * nrow.pt
  ntables <- nrow(data.mats$colmat)
    
  #
  # Pre-computation:
  #

  if (ncol(exitpoll)!=(ncol.pt*nrow.pt))
    stop("exitpoll must have nrow.pt*ncol.pt columns")
  if (nrow(exitpoll)!=ntables)
    stop("exitpoll must have the same number of rows as the data matrix")

  # Check there are at least 2 rows/cols:
  if (nrow.pt<2 || ncol.pt<2)
    stop("Table must be at least 2x2")

  # Check the NNtots is valid:
  if (any(apply(NNtots[,1:nrow.pt],1,sum) != apply(NNtots[,(nrow.pt+1):(nrow.pt+ncol.pt)],1,sum)))
    stop("Row and column totals do not sum to the same total total")

  # Handle the eschew argument:
  # We assume that if eschew=TRUE, then all of the matrices passed as arguments
  # including ncol.pt, already include the extra eschew column.
  if (eschew){

    # Error check:
    if (ncol.pt<3)
      stop("At least 3 columns required for eschew=TRUE")

    # Any other error checks?

  }

  mu.len <- nrow.pt*(ncol.pt-1)
  if (length(mu.vec.0)!=mu.len)
    stop("mu.vec.0 must be of length nrow.pt*(ncol.pt-1)")

	# Error check mu.vec.cu:
  if (length(mu.vec.cu)!=mu.len)
    stop("mu.vec.cu be of length nrow.pt*(ncol.pt-1)")

  # Can just use names from the standard formula:
  # This produces a vector of c("KKrn1cn1","KKrn1cn2",...,"KKrnRcnC"),
  # where rni and cnj are the names of row i and col j respectively.
  KKnames <- apply(expand.grid(colnames.local,rownames.local),1,
                   function(x){sprintf("KK%s%s",x[2],x[1])})

  # Or, can use dot-separated version:
  KKdotnames <- apply(expand.grid(colnames.local,rownames.local),1,
                      function(x){sprintf("KK.%s.%s",x[2],x[1])})

  # Must require the exitpoll to have column names,
  # if it didn't then have to assume they specified them
  # in the correct order i.e., c("KKrn1cn1","KKrn1cn2",...,"KKrnRcnC")
  if (is.null(colnames(exitpoll))){
    colnames(exitpoll) <- KKnames
  }

  # Possibly re-order the columns of the exitpoll into the desired format:
  trynames  <- try({exitpoll <- exitpoll[,KKnames]},silent=TRUE)
  if (class(trynames)=="try-error"){
    # Non dot-separated colnames failed, try dot-separated:
    trydotnames <- try({exitpoll <- exitpoll[,KKdotnames]},silent=TRUE)
    if (class(trydotnames)=="try-error")
      stop("Invalid column names of exitpoll, must be of the form 'KKrnicnj', or 'KK.rni.cnj'")
  }

  # Exit poll should now be in row-wise (byrow=TRUE) format:
  # The C code MUST have it specified in this way!
  # i.e. (x_11,x_12,...x_1C,x_21,...,x_RC)

  KKcolsums <- matrix(NA,nrow=nrow(NNtots),ncol=ncol.pt)
  for (i in 1:ncol.pt)
    KKcolsums[,i] <- apply(exitpoll[,(i+((0:(nrow.pt-1))*ncol.pt))],1,sum) 

  KKrowsums <- matrix(NA,nrow=nrow(NNtots),ncol=nrow.pt)
  for (i in 1:nrow.pt)
    KKrowsums[,i] <- apply(exitpoll[,(1+(i-1)*ncol.pt):(i*ncol.pt)],1,sum) 

  # NNtots has rowsums first, then colsums. We compare against this:
  # NOTE: This version of KKtots is not used beyond this error check.
  KKtots <- cbind(KKrowsums,KKcolsums)

  if (any(KKtots>NNtots)){
    cat("Some of the exit poll row/column totals exceeds the table row/col totals,\n")
    cat("Problems occur in the following precincts:\n")
    problem.rows <- c(1:nrow(KKtots))[apply(KKtots-NNtots,1,function(x){any(x>0)})]
    print(rownames(KKtots)[problem.rows])
    cat("The exit poll totals in these precincts are:\n")
    KKtotsprint <- KKtots[problem.rows,]
    rownames(KKtotsprint) <- rownames(KKtots)[problem.rows]
    cat("The count totals in these precincts are:\n")
    NNtotsprint <- NNtots[problem.rows,]
    rownames(NNtotsprint) <- rownames(NNtots)[problem.rows]
    stop("Invalid exit poll")
  }

  # Compute the bounds:   
  MMbounds <- gq.bounds(nrow.pt,ncol.pt,NNtots-KKtots)
  NNbounds <- MMbounds + cbind(exitpoll,exitpoll)

  # Force correct storage types:
  storage.mode(NNtots) <- "double"
  storage.mode(NNbounds) <- "double"
  storage.mode(MMbounds) <- "double"
  storage.mode(exitpoll) <- "double"
  storage.mode(mu.vec.cu) <- "double"
  NNtots    <- as.matrix(NNtots)
  NNbounds  <- as.matrix(NNbounds)
  MMbounds  <- as.matrix(MMbounds)
  exitpoll  <- as.matrix(exitpoll)
  mu.vec.cu <- as.matrix(mu.vec.cu)

  # Function to determine the choice of sr.reps:
  # x<0.11 => nchg
  # 0.11<=x <0.20 => 1 rep multinomial
  # 0.20<=x=<0.80 => 2 rep multinomial
  # x>0.80 => nchg
  #

  ep.rep.func <- function(x){
    r <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    for (i in 1:nrow(x)){
      for (j in 1:ncol(x)){
        tmp <- ifelse(x[i,j]<0.20,2,1)
        tmp <- ifelse((x[i,j]<0.11)||(x[i,j]>0.80),3,tmp)
        r[i,j] <- tmp
      }
    }
    return(r)
  }

  if (any(is.null(sr.probs))){
		# Default according to row fractions:
		rfracs <- NNtots[,1:nrow.pt]
		rfracs <- rfracs/matrix(rep(apply(rfracs,1,sum),nrow.pt),ncol=nrow.pt)
		sr.probs <- rfracs
	}
  storage.mode(sr.probs) <- "double"
  if (any(is.null(sr.reps))){
    # Default according to rep.func:
    sr.reps <- ep.rep.func(sr.probs)
  }
  storage.mode(sr.reps) <- "integer"
  
  # Is there a pre-specified starting state of NNs/MMs or random start?
  if (!is.null(NNs.start) || !is.null(MMs.start)){
				
    # User-supplied starting states. Error check them:
    storage.mode(NNs.start) <- "double"
    storage.mode(MMs.start) <- "double"
    NNs.start <- as.matrix(NNs.start)
    MMs.start <- as.matrix(MMs.start)

    # Error check the number of rows:
    if (!(nrow(NNs.start)==nrow(NNtots)))
      stop("Invalid number of rows in NNs.start")
    if (!(nrow(MMs.start)==nrow(NNtots)))
      stop("Invalid number of rows in MMs.start")

    # Error check the number of columns:
    if (!(ncol(NNs.start)==(nrow.pt*ncol.pt)))
      stop("Invalid number of cols in NNs.start")
    if (!(ncol(MMs.start)==(nrow.pt*ncol.pt)))
      stop("Invalid number of cols in MMs.start")

    # Error check for negative entries:
    if (any(NNs.start<0.0))
      stop("Entries of NNs.start must be non-negative")
    if (any(MMs.start<0.0))
      stop("Entries of MMs.start must be non-negative")

    # Pre-specified starting state for NNs:
    rstart.NNs.MMs <- FALSE

  } else {
    # Random starting state for NNs and MMs:
    rstart.NNs.MMs <- TRUE
  }
    
  # Is there a pre-specified starting state of the THETAS or random start?
  if (!(is.null(THETAS.start))){

    # User-supplied starting values. Error check:
    storage.mode(THETAS.start) <- "double"
    THETAS.start <- as.matrix(THETAS.start)
    if (!(nrow(THETAS.start)==nrow(NNtots)))
      stop("Invalid number of rows in THETAS.start")
    if (!(ncol(THETAS.start)==(nrow.pt*ncol.pt)))
      stop("Invalid number of cols in THETAS.start")
 
    # Pre-specified starting state for NNs:
    rstart.THETAS <- FALSE

  } else {
    # Random starting state for NNs:
    rstart.THETAS <- TRUE
  }

  # Error-check Diri:
  Diri <- rep(as.numeric(Diri),nrow(NNtots))
  if (length(Diri)!=nrow(NNtots))
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

    args <- list("NNtots" = NNtots,
                 "NNbounds" = NNbounds,
                 "ExitPoll" = exitpoll,
                 "MMbounds" = MMbounds,
                 "numrows_pt" = as.integer(nrow.pt),
                 "numcols_pt" = as.integer(ncol.pt),
                 "eschew" = as.integer(eschew),
                 "sr_probs" = as.matrix(sr.probs),
                 "sr_reps" = as.matrix(sr.reps),
                 "num_iters" = as.integer(num.iters),
                 "psi_0" = as.numeric(psi),
                 "kappa_0" = as.numeric(kappa),
                 "nu_0" = as.numeric(nu),
                 "nolocalmode" = as.numeric(nolocalmode),
                 "dof" = as.numeric(dof),
                 "mu_vec_0" = as.numeric(mu.vec.0),
                 "mu_vec_cu" = mu.vec.cu,
                 "use_Diri_every_vec" = as.numeric(Diri),
                 "rho_vec" = as.numeric(rho.vec),
                 "save_every" = as.numeric(save.every),
                 "numscans" = as.integer(numscans),
                 "rstart_NNs_MMs" = as.integer(rstart.NNs.MMs),
                 "NNs_start" = NNs.start,
                 "MMs_start" = MMs.start,
                 "rstart_THETAS" = as.integer(rstart.THETAS),
                 "THETAS_start" = THETAS.start,
                 "keep_restart_info" = as.integer(keep.restart.info),
                 "keepNNinternals" = as.double(keepNNinternals),
                 "keepTHETAS" = as.double(keepTHETAS),
                 "print_every" = as.integer(print.every),
                 "dbg" = as.integer(debug))
    gq.output <- .Call("AnalyzeWithExitPoll", PACKAGE="RxCEcolInf", args)

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
