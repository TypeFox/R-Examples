"wccsom" <- function(data, grid=somgrid(), rlen = 100,
                     alpha = c(0.05, 0.01),
                     radius = quantile(nhbrdist, 0.7),
                     init, nhbrdist, trwidth = 20, 
                     toroidal = FALSE, FineTune = TRUE,
                     keep.data = TRUE) 
{
  data <- as.matrix(data)
  datadims <- dimnames(data)
  dimnames(data) <- NULL
  nobj <- nrow(data)
  nvar <- ncol(data)
  nunits <- nrow(grid$pts)
    
  if (trwidth > 0)
    wghts <- 1 - (0:trwidth)/trwidth
  else
    wghts <- 1

  ## Default initialisation: random linear combinations
  ## of data points   
  if(missing(init)) {
    initwghts <- matrix(runif(nunits*nobj), nunits, nobj)
    initwghts <- sweep(initwghts, 1, rowSums(initwghts), FUN="/")
    init <- initwghts %*% data
  }
  
  ## In the C-code, it is easier to have data
  ## and codes transposed.
  data <- t(data)
  init <- t(init)
  codes <- init
  dimnames(codes) <- NULL

  ## declare space for arrays
  changes <- rep(0, rlen)
  acors  <- rep(0, nunits)
  data.acors <- rep(0, nobj)
  
  if (missing(nhbrdist))
    nhbrdist <- unit.distances(grid, toroidal)

  if (toroidal) {
    if (grid$topo == "hexagonal" & (grid$ydim %% 2 == 1))
      stop("Error: uneven number of rows (y) in hexagonal toroidal grid")
    radius <- radius*0.5
  }

  res <- .C("WCC_onlineSOM",
            data = as.double(data),
            codes = as.double(codes),
            nhbrdist = as.double(nhbrdist),
            alpha = as.double(alpha),
            radius = as.double(radius),
            trwdth = as.integer(trwidth),
            wghts = as.double(wghts),
            data.acors = as.double(data.acors),
            acors = as.double(acors),
            changes = as.double(changes),
            nobj = as.integer(nobj),
            nvar = as.integer(nvar),
            ncodes = as.integer(nunits),
            rlen = as.integer(rlen),
            PACKAGE = "wccsom")

  changes <- res$changes
  codes <- matrix(res$codes, nvar, nunits)
  acors <- res$acors
  data.acors <- res$data.acors
  
  if (FineTune) {
    res <- wcc.kmeans(data, data.acors, codes, acors, trwidth)
    acors <- res$acors
    codes <- res$codes
    if (keep.data) {
      unit.classif <- res$classif
      wccs <- res$wccs
    }
  } else {
    if (keep.data) {
      classif <- wccs <- rep(0, nobj)
      res <-  .C("wccassign",
                 data = as.double(data),
                 data.acors = as.double(data.acors),
                 codes = as.double(codes),
                 acors = as.double(acors),
                 trwdth = as.integer(trwidth),
                 wghts = as.double(wghts),
                 classif = as.integer(classif),
                 wccs = as.double(wccs),
                 as.integer(nobj),
                 as.integer(nvar),
                 as.integer(nunits),
                 PACKAGE = "wccsom")
      unit.classif <- res$classif
      wccs <- res$wccs
    }
  }

  if (keep.data) {
    structure(list(data = t(data), grid = grid, changes = changes,
                   codes = t(codes), trwdth = trwidth,
                   unit.classif = unit.classif, wccs = wccs,
                   data.acors = data.acors, acors = acors,
                   toroidal = toroidal, FineTune = FineTune),
              class = c("wccsom","SOM"))
  } else {
    structure(list(grid = grid, changes = changes, codes = t(codes),
                   trwdth = trwidth, acors = acors, toroidal = toroidal,
                   FineTune = FineTune), 
              class = c("wccsom","SOM"))
  }
}


### Calculate distances in a toroidal Kohonen map. Crude and
### slow implementation, but hey. When training a map takes 20 minutes,
### who cares about a few seconds?

unit.distances <- function(grid, toroidal)
{
  if (!toroidal) return(as.matrix(dist(grid$pts)))

  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[,2])/2
  
  result <- matrix(0, np, np)
  for (i in 1:(np-1)) {
    for (j in (i+1):np) {
      diffs <- abs(grid$pts[j,] - grid$pts[i,])
      if (diffs[1] > maxdiffx)
        diffs[1] <- 2*maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy)
        diffs[2] <- 2*maxdiffy - diffs[2]

      result[i,j] <- sum(diffs^2)
    }
  }

  sqrt(result + t(result))
}


## wcc.kmeans performs kmeans clustering using the WCC similarity
## criterion until convergence is reached
## data contains one object per column; data.acors is a vector of data
## autocorrelations, codes is a matrix of code vectors, and acors is a
## vector of code autocorrelations.

wcc.kmeans <- function(data, data.acors, codes, acors, tw, maxit=20) {
  wghts <- 1 - (0:tw)/tw
  nvar <- nrow(data) # data is transposed...
  nobj <- ncol(data)
  nunits <- length(acors)

  wccs <- oldclassif <- newclassif <- rep(0, nobj)
  if (missing(data.acors))
    data.acors <- wacmat(data, tw, wghts, do.transpose=FALSE)

  res <- .C("wccassign",
            data = as.double(data),
            data.acors = as.double(data.acors),
            codes = as.double(codes),
            acors = as.double(acors),
            as.integer(tw),
            wghts = as.double(wghts),
            oldclassif = as.integer(oldclassif),
            wccs = as.double(wccs),
            as.integer(nobj),
            as.integer(nvar),
            as.integer(nunits),
            PACKAGE = "wccsom")

  oldclassif <- res$oldclassif

  for (i in 1:maxit) {
    for (j in 1:nunits) {
### update codes: look at batchSom for another
###      possible implementation. Maybe faster.  
      whichones <- which(oldclassif == j)
      if (length(whichones)>0) {
        codes[,j] <- apply(data[,whichones, drop=FALSE], 1, mean)
      }
    }
    acors <- wacmat(codes, tw, wghts, do.transpose=FALSE)

    res <- .C("wccassign",
              data = as.double(data),
              data.acors = as.double(data.acors),
              codes = as.double(codes),
              acors = as.double(acors),
              as.integer(tw),
              wghts = as.double(wghts),
              newclassif = as.integer(newclassif),
              wccs = as.double(wccs),
              as.integer(nobj),
              as.integer(nvar),
              as.integer(nunits),
              PACKAGE = "wccsom")
    newclassif <- res$newclassif

    cat("\nIteration ", i, ": agreement ", sum(newclassif == oldclassif),
        " out of ", nobj, " cases", sep="")
    if (all(newclassif == oldclassif)) break

    oldclassif <- newclassif
  }

  cat("\n")
  
  if (i == maxit)
    warning("Maximum number of iterations reached in wcc.kmeans\n")
  
  list(codes = codes, acors = res$acors,
       classif = newclassif, wccs = res$wccs, i = i)
}



## Weighted Autocorrelation and Crosscorrelation functions

wac <- function(pattern1, trwdth, wghts)
{
  np <- length(pattern1)
  if (missing(wghts)) {
    if (trwdth > 0)
      wghts <- 1 - (0:trwdth)/trwdth
    else
      wghts <- 1
  }

  .C("wacdist",
     as.double(pattern1),
     as.integer(np),
     as.double(wghts),
     as.integer(trwdth),
     wacval = double(1),
     PACKAGE = "wccsom")$wacval
}


wacmat <- function(patterns, trwdth, wghts, do.transpose=TRUE)
{
  if (do.transpose)
    patterns <- t(patterns)
  
  nobj <- dim(patterns)[2]
  np <- dim(patterns)[1]
  
  if (missing(wghts)) {
    if (trwdth > 0)
      wghts <- 1 - (0:trwdth)/trwdth
    else
      wghts <- 1
  }

  wacval <- rep(0, nobj)
  .C("wacdists",
     as.double(patterns),
     as.integer(nobj),
     as.integer(np),
     as.double(wghts),
     as.integer(trwdth),
     wacval = as.double(wacval),
     PACKAGE = "wccsom")$wacval
}


wcc <- function(pattern1, pattern2, trwdth, wghts, acors)
{
  np <- length(pattern1)
  if (length(pattern2) != np)
    stop("Patterns of different length")
  
  if (missing(wghts))  {
    if (trwdth > 0)
      wghts <- 1 - (0:trwdth)/trwdth
    else
      wghts <- 1
  }
  
  if (missing(acors))
    acors <- c(wac(pattern1, trwdth, wghts),
               wac(pattern2, trwdth, wghts))

  res <- .C("wccdist",
            as.double(pattern1),
            as.double(pattern2),
            as.integer(np),
            as.double(wghts),
            as.integer(trwdth),
            crossterm = double(1),
            PACKAGE = "wccsom")$crossterm

  res/(acors[1]*acors[2])
}


