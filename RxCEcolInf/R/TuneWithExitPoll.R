`TuneWithExitPoll` <-  
           function(fstring, exitpoll, data=NULL, num.runs=12, num.iters=10000,
                    rho.vec=rep(0.05, ntables),
                    kappa=10, nu=(mu.dim+6), psi=mu.dim, 
                    mu.vec.0=rep(log((0.45/(mu.dim-1))/0.55), mu.dim),
                    mu.vec.cu=runif(mu.dim, -3, 0),
                    nolocalmode=50,
                    sr.probs=NULL, sr.reps=NULL, 
                    numscans=1, Diri=100, dof=4, debug=1)
{
  data.mats <- get.mat(fstring, data)
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
  
  # Check there are at least 2 rows/cols:
  if (nrow.pt<2 || ncol.pt<2)
    stop("Table must be at least 2x2")

  # Check the NNtots is valid:
  if (any(apply(NNtots[,1:nrow.pt],1,sum) != apply(NNtots[,(nrow.pt+1):(nrow.pt+ncol.pt)],1,sum)))
    stop("Row and column totals do not sum to the same total total")

  # Check the dimension of the exit poll:
  if (ncol(exitpoll)!=(ncol.pt*nrow.pt))
    stop("exitpoll must have nrow.pt*ncol.pt columns")
  if (nrow(exitpoll)!=ntables)
    stop("exitpoll must have the same number of rows as the data matrix")

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
  trynames    <- try({exitpoll <- exitpoll[,KKnames]},silent=TRUE)
  if (class(trynames)=="try-error"){
    # Non dot-separated colnames failed, try dot-separated:
    trydotnames <- try({exitpoll <- exitpoll[,KKdotnames]},silent=TRUE)
    if (class(trydotnames)=="try-error")
      stop("Invalid column names of exitpoll, must be of the form 'KKrnicnj', or 'KK.rni.cnj'")
  }

  # Exit poll should be in row-wise (byrow=TRUE) format:
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
    if (is.null(rownames(NNtots))){
      print(problem.rows)
    } else {
      print(rownames(NNtots)[problem.rows])
    }
    cat("The exit poll totals in these precincts are:\n")
    KKtotsprint <- KKtots[problem.rows,]
    rownames(KKtotsprint) <- rownames(KKtots)[problem.rows]
    print(KKtotsprint)
    cat("The count totals in these precincts are:\n")
    NNtotsprint <- NNtots[problem.rows,]
    rownames(NNtotsprint) <- rownames(NNtots)[problem.rows]
    print(NNtotsprint)
    stop("Invalid exit poll")
  }

  MMbounds <- gq.bounds(nrow.pt,ncol.pt,NNtots-KKtots)
  NNbounds <- MMbounds + cbind(exitpoll,exitpoll)

  if (any(is.null(sr.probs)))
    sr.probs <- matrix(1/nrow.pt,nrow=nrow(NNtots),ncol=nrow.pt)
  if (any(is.null(sr.reps)))
    sr.reps <- matrix(1,nrow=nrow(NNtots),ncol=nrow.pt)
  
  # Force correct storage types:
  storage.mode(NNtots) <- "double"
  storage.mode(NNbounds) <- "double"
  storage.mode(MMbounds) <- "double"
  storage.mode(exitpoll) <- "double"
  storage.mode(mu.vec.cu) <- "double"
  storage.mode(sr.probs) <- "double"
  storage.mode(sr.reps) <- "integer"
  NNtots    <- as.matrix(NNtots)
  NNbounds  <- as.matrix(NNbounds)
  MMbounds  <- as.matrix(MMbounds)
  exitpoll  <- as.matrix(exitpoll)
  mu.vec.cu <- as.matrix(mu.vec.cu)
  sr.probs  <- as.matrix(sr.probs)
  sr.reps   <- as.matrix(sr.reps)

  # Error-check Diri:
  Diri <- rep(as.numeric(Diri),nrow(NNtots))
  if (length(Diri)!=nrow(NNtots))
    stop("Diri should be a scalar quantity")
  if (any(Diri<0.0))
    stop("Diri should be non-negative")
  storage.mode(Diri) <- "double"

  args <- list("NNtots" = NNtots,
               "NNbounds" = NNbounds,
               "ExitPoll" = exitpoll,
               "MMbounds" = MMbounds,
               "numrows_pt" = as.integer(nrow.pt),
               "numcols_pt" = as.integer(ncol.pt),
               "sr_probs" = sr.probs,
               "sr_reps" = sr.reps,
               "num_iters" = as.integer(num.iters),
               "psi_0" = as.numeric(psi),
               "kappa_0" = as.numeric(kappa),
               "nu_0" = as.numeric(nu),
               "nolocalmode" = as.numeric(nolocalmode),
               "num_runs" = as.integer(num.runs),
               "dof" = as.numeric(dof),
               "num_runs" = as.integer(num.runs),
               "mu_vec_0" = as.numeric(mu.vec.0),
               "mu_vec_cu" = as.numeric(mu.vec.cu),
               "use_Diri_every_vec" = as.numeric(Diri),
               "rho_vec" = as.numeric(rho.vec),
               "numscans" = as.integer(numscans),
               "dbg" = as.integer(debug))
  .Call("TuneWithExitPoll", PACKAGE="RxCEcolInf", args)
}
