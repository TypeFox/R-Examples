### simSummary.R
###-----------------------------------------------------------------------------

simSummary <- structure(

  ## --- Function ---

function # Simulation summary
  ### \code{simSummary} eases the process of summarizing simulation results.
  ### Simulations often produce some intermediate results (some focal statistic(s)),
  ### that need to be summarized over many simulation replicates. \code{simSummary}
  ### helps with summarizing these focal statistics.
(
  x,   ##<< an (outer) list of (inner) lists, where each inner list has exactly the same structure (see examples)
  FUN=c("length", "nobs", "mean", "sd", "min", "max"), ##<< character, summary statistics function names
  ...  ##<< arguments passed to summary functions
) {
  
  ##details<< \code{simSummary} accepts as input an (outer) list of (inner) lists,
  ## where all inner lists must have the same structure and only scalars, vectors,
  ## matrices, and arrays can be used in inner lists. Function combines all inputs
  ## in a list of arrays and summarizes array values with specified functions that
  ## can work on vector like inputs.

  ## --- Constants ---

  nS <- length(x)      ## number of simulations
  nO <- length(x[[1]]) ## number of objects in each simulation
  nF <- length(FUN)    ## number of summarizing functions

  ## --- Test input (x) ---

  tmp <- "'x' must be a list"
  if(!is.list(x)) {
    stop(tmp)
  } else {
    if(is.data.frame(x)) {
      stop(tmp)
    }
  }

  test <- sapply(X=x, FUN=is.list) 
  if(any(!test)) stop("all elements of 'x' must be a list")

  test <- sapply(X=x[[1]], FUN=function(z) is.numeric(z) & (is.vector(z) | is.matrix(z) | is.array(z)))
  if(any(!test)) {
    stop(paste("elements of inner lists must be either numeric vector, matrix, or array: ", paste(names(x[[1]])[!test], collapse=", "), sep=""))
  }

  dimFun <- function(z)
  {
    ret <- dim(z)
    if(is.null(ret)) {
      ret <- length(z)
    }
    ret
  }
  testDim <- lapply(X=x, FUN=function(z) {lapply(X=z, FUN=dimFun)})
  for(i in 2:nS) {
    ## i <- 2
    test <- !identical(testDim[[i]], testDim[[i-1]])
    if(test) stop("all elements of 'x' must have exactly the same structure")
  }

  myDim <- testDim[[1]]
  nD <- sapply(testDim[[1]], length) ## dimension lengths

  ## --- Test input (FUN) ---

  tmp <- 1:10
  ret <- list()
  for(f in 1:nF) {
    ## f <- 1
    ret[[f]] <- do.call(what=FUN[f], args=list(tmp))
  }
  test <- sapply(X=ret, FUN=function(z) is.vector(z) & length(z) == 1)
  if(any(!test)) {
    stop(paste("summary functions must return a single (scalar) value: ", paste(FUN[!test], collapse=", "), sep=""))
  }

  ## --- Setup working objects --- 

  arr <- vector(mode="list", length=nO)
  names(arr) <- names(x[[1]])
  ret <- vector(mode="list", length=nF)
  names(ret) <- FUN
  for(i in 1:nF) {
    ret[[i]] <- x[[1]]
  }

  ## --- Build arrays ---

  if(FALSE) {

    ## Structure
    for(o in 1:nO) {
      ## o <- 1
      arr[[o]] <- array(data=NA, dim=c(testDim[[1]][[o]], nS))
    }

    ## Content
    for(s in 1:nS) {
      ## s <- 1
      for(o in 1:nO) {
        ## o <- 1
        tmp <- dim(arr[[o]])
        tmp[length(tmp)] <- s
        ## afill(arr[[o]][s], ?, ?) <- x[[s]][[o]]
      }
    }
  }

  ## Content
  ## - NOTE: the cost of this approach is nO + nS*nO
  ## - BEWARE: this approach uses reallocation for each abind - might be slow!
  for(o in 1:nO) {
    arr[[o]] <- x[[1]][[o]]
  }
  for(s in 2:nS) {
    ## s <- 2
    for(o in 1:nO) {
      ## o <- 1
      arr[[o]] <- abind(arr[[o]], x[[s]][[o]], along=nD[o]+1)
    }
  }

  ## --- Summarize ---

  for(i in 1:nO) {
    myDim[[i]] <- 1:nD[[i]]
  }

  ## - NOTE: the cost of this approach is nF*nO*?
  for(f in 1:nF) {
    ## f <- 1
    for(o in 1:nO) {
      ## o <- 1
      ret[[f]][[o]] <- apply(X=arr[[o]], MARGIN=myDim[[o]], FUN=FUN[f], ...)
    }
  }

  ## --- Return ---

  ret

  ##value<< The return element is also a list (outer) of lists (inner), where
  ## each inner list has the same structure as inner lists of input, but holding
  ## one of the summary statistics - one summary statistics per one inner list.

},## --- Function end ---

  ## --- Examples ---

ex=function() {

  ## Create simple input from a rather silly simulation
  simFun <- function(x)
  {
    ret <- list()
    ret$s <- rnorm(n=1)
    ret$v <- rnorm(n=5)
    ret$m <- matrix(rnorm(n=5*5), nrow=5, ncol=5)
    ret$a <- array(rnorm(n=4*3*2), dim=c(4, 3, 2))
    ret
  }
  sim <- list()
  sim$sim1 <- simFun()
  sim$sim2 <- simFun(x=0)
  sim$sim3 <- simFun(x=1)

  ## Simulation summary (just mean and standard deviation)
  simSummary(x=sim, FUN=c("mean", "sd"))

  ## Can handle simulations in process too = handle NA values
  sim$sim3$s <- NA
  sim$sim3$v[5] <- NA
  simSummary(x=sim, FUN="mean")
  simSummary(x=sim, FUN="mean", na.rm=TRUE)

  ## Unit tests (automatic run elsewhere)
  ## summary(runTest(test(simSummary)))

})## --- Examples end ---

  ## --- Unit tests ---

test(simSummary) <- function()
{

  ## Create simple input from a rather silly simulation - fixed to ease testing
  simFun <- function(x)
  {
    ret <- list()
    ret$s <- 1 + x
    ret$v <- 1:5 + x
    ret$m <- matrix(1:25, nrow=5, ncol=5) + x
    ret$a <- array(1:(4*3*2), dim=c(4, 3, 2)) + x
    ret
  }
  sim <- list()
  sim$sim1 <- simFun(x=-1)
  sim$sim2 <- simFun(x=0)
  sim$sim3 <- simFun(x=1)

  ## ... with NA
  simNA <- sim
  simNA$sim3$s          <- NA
  simNA$sim3$v[5]       <- NA
  simNA$sim3$m[5, 5]    <- NA
  simNA$sim3$a[4, 3, 2] <- NA

  ## Test input - 'x' must be a list
  (checkException(simSummary(x=1)))
  (checkException(simSummary(x=data.frame(1:3, 1:3))))

  ## Test input - "all elements of 'x' must be a list"
  tmp <- sim; tmp$sim4 <- 1
  (checkException(simSummary(x=tmp)))

  ## Test input - "elements of inner lists must be either numeric vector, matrix, or array"
  tmp <- sim; tmp$sim1$l <- list(1:10)
  (checkException(simSummary(x=tmp)))
  tmp <- sim; tmp$sim1$d <- data.frame(a=1:10, b=1:10)
  (checkException(simSummary(x=tmp)))
  tmp <- sim; tmp$sim1$s <- "a"
  (checkException(simSummary(x=tmp)))

  ## Test input - all elements of 'x' must have exactly the same structure
  tmp <- sim; tmp$sim1$s <- 1:3
  (checkException(simSummary(x=tmp)))

  ## Test input - "summary functions must return a single (scalar) value"
  (checkException(simSummary(x=tmp, FUN=range)))

  ## Test summary values - length
  ret <- simSummary(x=sim, FUN="length")
  (checkEquals(ret$length$s, 3, checkNames=FALSE))
  (checkEquals(ret$length$v, rep(3, 5), checkNames=FALSE))
  tmp <- ret$length$m; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$length$m, tmp, checkNames=FALSE))
  tmp <- ret$length$a; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$length$a, tmp, checkNames=FALSE))
  
  ## Test summary values - nobs
  ret <- simSummary(x=sim, FUN="nobs")
  (checkEquals(ret$nobs$s, 3, checkNames=FALSE))
  (checkEquals(ret$nobs$v, rep(3, 5), checkNames=FALSE))
  tmp <- ret$nobs$m; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$nobs$m, tmp, checkNames=FALSE))
  tmp <- ret$nobs$a; tmp[] <- 3 # to avoid attributes issues
  (checkEquals(ret$nobs$a, tmp, checkNames=FALSE))

  ## Test summary values - nobs (in real action)
  ret <- simSummary(x=simNA, FUN="nobs")
  (checkEquals(ret$nobs$s, 2, checkNames=FALSE))
  (checkEquals(ret$nobs$v, c(rep(3, 4), 2), checkNames=FALSE))
  tmp <- ret$nobs$m; tmp[] <- 3; tmp[5, 5] <- 2 # to avoid attributes issues
  (checkEquals(ret$nobs$m, tmp, checkNames=FALSE))
  tmp <- ret$nobs$a; tmp[] <- 3; tmp[4, 3, 2] <- 2 # to avoid attributes issues
  (checkEquals(ret$nobs$a, tmp, checkNames=FALSE))

  ## Test summary values - mean
  ret <- simSummary(x=sim, FUN="mean")
  (checkEquals(ret$mean$s, 1, checkNames=FALSE))
  (checkEquals(ret$mean$v, 1:5, checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- 1:25 # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- 1:24 # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - mean & NA
  ret <- simSummary(x=simNA, FUN="mean")
  (checkEquals(ret$mean$s, as.numeric(NA), checkNames=FALSE))
  (checkEquals(ret$mean$v, c(1:4, NA), checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- c(1:24, NA) # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- c(1:23, NA) # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - mean & NA but with na.rm=TRUE
  ret <- simSummary(x=simNA, FUN="mean", na.rm=TRUE)
  (checkEquals(ret$mean$s, 0.5, checkNames=FALSE))
  (checkEquals(ret$mean$v, c(1:4, 4.5), checkNames=FALSE))
  tmp <- ret$mean$m; tmp[] <- c(1:24, 24.5) # to avoid attributes issues
  (checkEquals(ret$mean$m, tmp, checkNames=FALSE))
  tmp <- ret$mean$a; tmp[] <- c(1:23, 23.5) # to avoid attributes issues
  (checkEquals(ret$mean$a, tmp, checkNames=FALSE))

  ## Test summary values - min
  ret <- simSummary(x=sim, FUN="min")
  (checkEquals(ret$min$s, 0, checkNames=FALSE))
  (checkEquals(ret$min$v, 0:4, checkNames=FALSE))
  tmp <- ret$min$m; tmp[] <- 0:24 # to avoid attributes issues
  (checkEquals(ret$min$m, tmp, checkNames=FALSE))
  tmp <- ret$min$a; tmp[] <- 0:23 # to avoid attributes issues
  (checkEquals(ret$min$a, tmp, checkNames=FALSE))

  ## Test summary values - max
  ret <- simSummary(x=sim, FUN="max")
  (checkEquals(ret$max$s, 2, checkNames=FALSE))
  (checkEquals(ret$max$v, 2:6, checkNames=FALSE))
  tmp <- ret$max$m; tmp[] <- 2:26 # to avoid attributes issues
  (checkEquals(ret$max$m, tmp, checkNames=FALSE))
  tmp <- ret$max$a; tmp[] <- 2:25 # to avoid attributes issues
  (checkEquals(ret$max$a, tmp, checkNames=FALSE))

} ## --- Unit tests end ---

###-----------------------------------------------------------------------------
### simSummary.R ends here
