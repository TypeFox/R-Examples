if(0){ # Delete this whole file (gendataExitPoll.R) when Cassandra
  ###     gets done with her validation simulations,
  ###     use the gendata.ep function in SimDataEP.R instead


##  Function to generate test data from the model, include exit poll

gendata <- function(nprecincts = 175, nrowcat = 3, ncolcat = 3,
                    rowcatnames = c("bla", "whi", "his"),
                    colcatnames = c("Dem", "Rep", "Abs"),
                    nprecincts.ep = 0,
                    samplefrac.ep = .1,
                    mu0 = c(0, -2, -2, 0, -.5, -1.5),
                    K0 = diag(.3, 6),
                    nu0 = 12, Psi0 = diag(6)*6,
                    lambda = 1000,
                    alpha = c(.4, .5, .2) * 2,
                    outfile=NULL){


  
  ## #####################################################
  ## utility functions
  
  "rdirichlet" <- function (n, alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
  }
  
  "rwish" <-
    function(v, S) {
      if (!is.matrix(S))
        S <- matrix(S)
      if (nrow(S) != ncol(S)) {
        stop(message="S not square in rwish().\n")
      }
      if (v < nrow(S)) {
        stop(message="v is less than the dimension of S in rwish().\n")
      }
      p <- nrow(S)
      CC <- chol(S) 
      Z <- matrix(0, p, p)
      diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
      if(p > 1) {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
      }
      return(crossprod(Z %*% CC))
    }
  
  
  "riwish" <-
    function(v, S) {
      return(solve(rwish(v,solve(S))))
    }
  ## #####################################################




  
  
  ## #####################################################
  ## create storage matrices
  ## mydata is marginal total data
  mydata <- matrix(NA, nrowcat+ncolcat+1, nprecincts)
  rownames(mydata) <- c(rowcatnames, colcatnames, "nunits")
  colnames(mydata) <- paste("precinct", 1:nprecincts, sep="")

  ## omegamat is the matrix of omega_i values (i on the columns)
  omegamat <- matrix(NA, (ncolcat-1)*nrowcat, nprecincts)
  rownames(omegamat) <- paste(
                              rep(
                                  paste(colcatnames[1:(ncolcat-1)],
                                        colcatnames[ncolcat], sep="."),
                                  nrowcat),
                              rep(rowcatnames, rep(ncolcat-1, nrowcat)),
                              sep=".")
  
  
  ## completedata.list is a list of the complete data table
  completedata.list <- vector("list", nprecincts)
  ## #####################################################


  
  


  
  ## generate hyper-parameters (Level 3)
  mu <- rmvnorm(1, mean=mu0, sigma=K0)
  Sigma <- riwish(nu0, Psi0) ## may need to be changed depending on
                             ## parameterization of inv. Wishart

  ## add row and column names to parameters
  rownames(K0) <- colnames(K0) <- rownames(omegamat)
  colnames(mu) <- rownames(omegamat)
  rownames(Sigma) <- colnames(Sigma) <- rownames(omegamat)

  for (i in 1:nprecincts){

    ## generate Level 2 parameters
    omega.i <- rmvnorm(1, mean=mu, sigma=Sigma)


    ## convert omega.i to probabilities (theta.i)
    omega.i.mat <- matrix(omega.i, ncolcat-1, nrowcat)
    omega.i.aug.mat <- rbind(omega.i.mat, 0)
    theta.i.mat <- matrix(NA, nrow(omega.i.aug.mat), ncol(omega.i.aug.mat))
    for (j in 1:ncol(theta.i.mat)){
      theta.i.mat[,j] <- exp(omega.i.aug.mat[,j]) /
        sum( exp( omega.i.aug.mat[,j]))
    }

    
    ## generate number of units in table i
    nunits.i <- rpois(1, lambda)

    ## generate the row totals
    rowtotal.probs <- rdirichlet(1, alpha)
    rowtotals <- rmultinom(1, nunits.i, rowtotal.probs)

    ## generate the interior of the table
    mytable <- matrix(NA, nrowcat, ncolcat)
    rownames(mytable) <- rowcatnames
    colnames(mytable) <- colcatnames
    for (j in 1:nrowcat){
      mytable[j,] <- rmultinom(1, rowtotals[j], theta.i.mat[,j])
    }
    coltotals <- apply(mytable, 2, sum)

    
    ## fill out storage matrices
    mydata[,i] <- c(rowtotals, coltotals, sum(rowtotals))
    omegamat[,i] <- omega.i
    completedata.list[[i]] <- mytable
    
  } ## end i in 1:nprecincts

  

  ##  create exit poll simulation if desired
    #  set up the exit poll matrix
    #     exit poll matrix has one row for each precinct
    #     row is all 0s if precinct was not in exit poll sample
    #     format is r1c1, r1c2, . . . r1cC, r2c2, . . . . . . rRcC
  returnmat.ep <- matrix(0, nrow = nprecincts,
                         ncol = nrowcat * ncolcat)
  if(nprecincts.ep){

    #  make the columnames
    colnames.ep <- rep(NA, ncol(returnmat.ep))
    counter <- 1
    for(ii in 1:nrowcat){
      for(jj in 1:ncolcat){
        colnames.ep[counter] <- paste("KK.", rowcatnames[ii], ".",
                                      colcatnames[jj], sep = "")
        counter <- counter + 1
      }
    }
    colnames(returnmat.ep) <- colnames.ep

    #  get precincts in sample, srs of precincts
    if(nprecincts.ep > nprecincts) stop("Number of exit poll precincts greater than number of precincts.")
    sampprecincts.ep <- sample((1:nprecincts),
                                 nprecincts.ep, replace = F)

    #  obtain the samples, not very elegant but should work so long
    #      as samplefrac.ep isn't too big
    for(ii in 1:nprecincts.ep){
      precnum <- sampprecincts.ep[ii]
      truecountvec <- as.vector(t(completedata.list[[precnum]]))
      prectot <- sum(truecountvec)
      nnsamp <- floor(prectot * samplefrac.ep)
      draw <- rep(prectot, length(truecountvec))
      while(min(truecountvec - draw) < 0){
        draw<- rmultinom(1, nnsamp, truecountvec/nnsamp)
      }
      returnmat.ep[precnum,] <- draw
    }
  }

  ##  put mydata in format corresponding to new Analyze formula interface
  mydata <- t(mydata)
  tempcolct <- ncol(mydata)
  mydata <- mydata[,1:(tempcolct-1)]
  
  ## create output object
  GQ.sim.object <- list(GQdata=mydata, GQexitpoll = returnmat.ep,
                        omega.matrix=omegamat,
                        interior.tables=completedata.list,
                        mu0=mu0, K0=K0, nu0=nu0, Psi0=Psi0,
                        mu=mu, Sigma=Sigma,
                        Sigma.diag=sqrtdiag(Sigma),
                        Sigma.cor=cov2cor(Sigma),
                        nprecincts=nprecincts,
                        nrowcat=nrowcat, ncolcat=ncolcat,
                        rowcatnames=rowcatnames, colcatnames=colcatnames,
                        lambda=lambda, alpha=alpha)
  
  if (!is.null(outfile)){
    save(GQ.sim.object, file=outfile)
  }
  
  return(GQ.sim.object)

}


}
#  order of output in current version of Analyze.R:
#     mus
#     sds of SIGMA
#     correlations of SIGMA, all first col, all second col, etc.
#     NNs, row major along precinct table
#     LAMBDAs, row major along precinct table
#     TURNOUT, in obvious order
#     GAMMAs, in obvious order
#     BETAs, row major along precinct table, but skipping last column
