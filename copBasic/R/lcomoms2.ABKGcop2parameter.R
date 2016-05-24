"lcomoms2.ABKGcop2parameter" <-
function(solutionenvir=NULL,
         T2.12=NULL, T2.21=NULL,
         T3.12=NULL, T3.21=NULL,
         T4.12=NULL, T4.21=NULL,
         t2eps=0.1, t3eps=0.1, t4eps=0.1,
         compt2erruv=TRUE, compt2errvu=TRUE,
         compt3erruv=TRUE, compt3errvu=TRUE,
         compt4erruv=TRUE, compt4errvu=TRUE,
         uset3err=TRUE, uset4err=FALSE,
         setreturn=1, maxtokeep=1e5) {

  if(is.null(T2.12)) {
    warning("T2.12 is NULL")
    return()
  }
  if(is.null(T2.21)) {
    warning("T3.21 is NULL")
    return()
  }
  if(is.null(T3.12)) {
    warning("T3.12 is NULL")
    return()
  }
  if(is.null(T3.21)) {
    warning("T3.21 is NULL")
    return()
  }
  if(is.null(T4.12)) {
    warning("T4.12 is NULL")
    return()
  }
  if(is.null(T4.21)) {
    warning("T4.21 is NULL")
    return()
  }

  if(is.null(solutionenvir)) {
    warning("Need the solution environment")
    return()
  }

  SOLUTIONS <- matrix(ncol=12, nrow=maxtokeep)

  eachfile <- ls(solutionenvir)
  eachfile <- sample(eachfile, length(eachfile), replace=FALSE)
  n <- 0; i <- 0

  # need the column names for a much later attributing
  tmp.lcom <- get(eachfile[1], envir=solutionenvir)
  if(is.data.frame(tmp.lcom)) {
    the.names <- names(tmp.lcom)
  } else {
    the.names <- unlist(attributes(tmp.lcom)$dimnames[2])
  }

  # For each "file" of solutions in the given environment
  # extract those arbitrary solutions meeting the *eps values
  # to build a much smaller suite of solutions from which we
  # will extract the solution with the smallest error
  for(file in eachfile) {
    if(i >= maxtokeep) next;
    message("Searching in ", file, appendLF=FALSE)
    # lcom contains the simulated L-comoments from known parameters
    lcom <- get(file, envir=solutionenvir)
    #print(lcom)
    if(is.data.frame(lcom)) lcom <- as.matrix(lcom)
    ix <- 1:attributes(lcom)$dim[1]
    ixn <- length(ix)

    n <- n + length(ix)


    # abs sum of the two offsets of 1 wrt 2 and 2 wrt 1
    ifelse(compt2erruv, t2erruv <- abs(T2.12 - lcom[,5]), t2erruv <- rep(0, ixn))
    ifelse(compt2errvu, t2errvu <- abs(T2.21 - lcom[,6]), t2errvu <- rep(0, ixn))
    t2err <- t2erruv + t2errvu
    t2errmin <- min(t2err) # capture the minimum for filtering

    # now create sub index of only those small errors in Tau2
    ixt2  <- ix[abs(t2err - t2errmin) < t2eps]
    ixt2n <- length(ixt2)

    # compute another abs sum of the two offsets of 1 wrt 2 and 2 wrt 1
    ifelse(compt3erruv, t3erruv <- abs(T3.12 - lcom[ixt2,7]), t3erruv <- rep(0, ixt2n))
    ifelse(compt3errvu, t3errvu <- abs(T3.21 - lcom[ixt2,8]), t3errvu <- rep(0, ixt2n))

    t3err <- t3erruv + t3errvu
    t3errmin <- min(t3err) # capture the minimum for filtering

    # now create sub index of only those small errors in Tau3
    ixt3  <- ixt2[abs(t3err - t3errmin) < t3eps]
    ixt3n <- length(ixt3)
    if(ixt3n == 0) {
       warning("No potential solutions in L-coskew, resetting counts to L-correlation")
       ixt3  <- ixt2
       ixt3n <- ixt2n
    }
    # compute another abs sum of the two offsets of 1 wrt 2 and 2 wrt 1
    ifelse(compt4erruv, t4erruv <- abs(T4.12 - lcom[ixt3,9]), t4erruv <- rep(0, ixt3n))
    ifelse(compt4errvu, t4errvu <- abs(T4.21 - lcom[ixt3,10]), t4errvu <- rep(0, ixt3n))
    t4err <- t4erruv + t4errvu
    t4errmin <- min(t4err) # capture the minimum for filtering

    # now create sub index of only those small errors in Tau4
    ixt4 <- ixt3[abs(t4err - t4errmin) < t4eps]
    ixt4n <- length(ixt4)
    if(ixt4n == 0) {
       warning("No potential solutions in L-cokurtosis, resetting counts to L-coskew")
       ixt4  <- ixt3
       ixt4n <- ixt3n
    }
    # Now pass through each of the potential solutions from this
    # file and store them in the SOLUTIONS matrix, which will be used
    # later to extract the minimum solution
    k <- length(lcom[ixt4,1])
    message(" found incremental solutions = ",k)

    if(i+k > maxtokeep) {
       k <- maxtokeep - i
    }
    for(j in 1:k) {
       if(i+j > maxtokeep) next;
       SOLUTIONS[i+j, ] <- sapply(1:12, function(x) { return(lcom[ixt4[j],x]) } )
    }
    i <- i + k # increment the file number
    #print(SOLUTIONS)
  }

  #message("Previous simulations checked: ",n)
  S <- SOLUTIONS[complete.cases(SOLUTIONS),] # safety check on NA, but
  #message("Truncating solutions in case of NAs")
  S <- as.matrix(S)
  #message("Solutions converted by as.matrix(S), just in case")
  if(attributes(S)$dim[2] == 1) S <- t(S)

  #print(attributes(S))
  # should have already been done by the preprocessing of the solutionenvir
  ns <- attributes(S)$dim[1] # number of incoming potential solutions
  nc <- attributes(S)$dim[2] # number of columns

  if(ns == 0) {
     warning("No solutions found")
  }

  message("Potential solutions found: ",ns)
  S <- as.data.frame(S)
  names(S) <- the.names
  S <- as.list(S)

  # Compute the 6 residuals
  S$t2.12res <- S$T2.12 - T2.12
  S$t2.21res <- S$T2.21 - T2.21
  S$t3.12res <- S$T3.12 - T3.12
  S$t3.21res <- S$T3.21 - T3.21
  S$t4.12res <- S$T4.12 - T4.12
  S$t4.21res <- S$T4.21 - T4.21

  # Rescale the abs of the residuals to put them all on the
  # same relative unit footing for combining
  ifelse(compt2erruv, t2erruv <- scale(abs(S$t2.12res), center=FALSE), t2erruv <- rep(0, ns))
  ifelse(compt2errvu, t2errvu <- scale(abs(S$t2.21res), center=FALSE), t2errvu <- rep(0, ns))
  S$t2err <- t2erruv + t2errvu

  ifelse(compt3erruv, t3erruv <- scale(abs(S$t3.12res), center=FALSE), t3erruv <- rep(0, ns))
  ifelse(compt3errvu, t3errvu <- scale(abs(S$t3.21res), center=FALSE), t3errvu <- rep(0, ns))
  S$t3err <- t3erruv + t3errvu

  ifelse(compt4erruv, t4erruv <- scale(abs(S$t4.12res), center=FALSE), t4erruv <- rep(0, ns))
  ifelse(compt4errvu, t4errvu <- scale(abs(S$t4.21res), center=FALSE), t4errvu <- rep(0, ns))
  S$t4err <- t4erruv + t4errvu

  S$jointerr <- S$t2err # L-correlation is always desired in minimization

  nb <- 2 # this is a counter used in a later sum() to return the
  # error of the fit back into the native units of the L-comoments
  # the value is 2 because of 1 wrt 2 and 2 wrt 1 with L-comoments
  if(uset3err) {
    nb <- nb + 2
    S$jointerr <- S$jointerr + S$t3err
  }
  if(uset4err) {
    nb <- nb + 2
    S$jointerr <- S$jointerr + S$t4err
  }

  S  <- as.data.frame(S)
  ix <- order(S$jointerr) # ordering indices of the solutions
  S  <- S[ix,] # order the entire data frame
  row.names(S) <- 1:length(S[,1]) # reset the row numbers
  i <- setreturn # this permits the user to "see"
  # other nearby solutions with the ith index used in the next line
  z <- as.list(S[i,1:nc]) # preparing the list object for the return()

  # Create a table for attachment to the return'ed list of real residuals
  # ** not rescaled units **
  resid <- data.frame(resT2.12=z$T2.12 - T2.12,
                      resT2.21=z$T2.12 - T2.12,
                      resT3.12=z$T3.12 - T3.12,
                      resT3.21=z$T3.21 - T3.21,
                      resT4.12=z$T4.12 - T4.12,
                      resT4.21=z$T4.21 - T4.21)
  my.names <- names(resid)
  resid[1,7] <- sum(abs(resid[1,1:nb])) # compute the true native unit
  # error of the L-comoments to the fitted copula
  names(resid) <- c(my.names, "SUM_ABS_USED_RESIDUALS")
  z$residuals <- resid

  z$solutions <- new.env() # the environment to house ALL of the
  # potential solutions meeting the *eps requirements
  assign("solutions", S, envir=z$solutions) # store the solutions

  return(z)
}
