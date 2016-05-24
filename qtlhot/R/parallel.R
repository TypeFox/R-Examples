## need to get rid of scandrop, zero, use highlod
######################################################################
# parallel.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: parallel.qtlhot, parallel.error, parallel.message,
#           qtlhot.phase0, qtlhot.phase1, qtlhot.phase2, qtlhot.phase3
######################################################################

######################################################################
## Directory path is indicated with "dirpath" character string
##
## I/O files are in internal R format (RData) using load/save commands
## except "groups.txt", which is used to set width of parallelization.
##
## I/O files are on AFS in
## /u/y/a/yandell/public/html/sysgen/qtlhot/condor
## or equivalently via URL
## http://www.stat.wisc.edu/~yandell/sysgen/qtlhot/condor
##
####################################################################################
parallel.qtlhot <- function(x, data = 1, ..., dirpath = ".")
{
  switch(x,
         qtlhot.phase1(dirpath, data, ...),
         qtlhot.phase2(dirpath, data, ...),
         qtlhot.phase3(dirpath, data, ...),
         parallel.error(1, x, data))
  
  parallel.error(0, x, data)
}
parallel.error <- function(num, phase = 0, index = 1)
{
  ## See file errorcodes.txt for explanation.
  if(phase == 5)
    outname <- "RESULT"
  else
    outname <- paste("RESULT", phase, index, sep = ".")
  write.table(num, file = outname,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  if(num)
    stop(parallel.message(num))

  invisible()
}
parallel.message <- function(num)
{
  msg <- as.matrix(c("OK",
                     "phase number not between 1 and 3",
                     "index to groups must be supplied as line number for file groups.txt",
                     "index to groups must be integer (line number for file groups.txt)",
                     "index to groups must be between 1 and number of lines in groups.txt",
                     "no bic RData files found",
                     "index to MCMC runs must be supplied",
                     "index to groups must be integer",
                     "index to groups must be between 1 and nruns",
                     "no MCMC RData files found",
                     "number of MCMC runs does not match nruns",
                     "number of BIC runs does not match groups"))
  dimnames(msg)[[1]] <- seq(0, nrow(msg) - 1)
  if(missing(num))
    write.table(msg, file = "errorcodes.txt",
                quote = FALSE, col.names = FALSE)
  else
    msg[as.character(num), 1]
}
####################################################################################
qtlhot.phase0 <- function(dirpath, init.seed = 92387475,
                          len = rep(400,16),
                          n.mar = 185, n.ind = 112,
                          n.phe = 100, latent.eff = 0, res.var = 1,
                          lod.thrs = c(4.63,4.17,3.93,3.76,3.65,3.56,3.47,3.39,3.34,3.30),
                          ...)
{
  ## PHASE 0: Cross object initialization. Create file "cross.RData" for Phase1,Phase3.
  
  set.seed(init.seed)

  mymap <- sim.map(len = len, n.mar = n.mar, include.x = FALSE, eq.spacing = TRUE)
  cross <- sim.cross(map = mymap, n.ind = n.ind, type = "bc")
  cross <- calc.genoprob(cross, step=0)
  
  ## Simulate new phenotypes.
  cross <- sim.null.pheno.data(cross, n.phe, latent.eff, res.var)
  
  save(cross, lod.thrs, file = file.path(dirpath, "cross.RData"), compress = TRUE)
}
####################################################################################
qtlhot.phase1 <- function(dirpath, index = 0,
                          params.file = "params.txt",
                          cross.file = "cross.RData",
                          cross.name = "cross",
                          n.ind = nind(cross),
                          lod.thrs = c(4.63,4.17,3.93,3.76,3.65,3.56,3.47,3.39,3.34,3.30),
                          alpha.levels = c(1:10) / 100, ## Nominal significance levels.
                          Nmax = 100, ## Maximum hotpsot size recorded.
                          n.perm = 1000, ## Number of permutations.
                          n.split = 100, ## Number of splits of permutations.
                          latent.eff = 1.5, ## Latent effect determines correlation among traits.
                          res.var = 1, ## Residual variance.
                          n.phe = nphe(cross), ## Number of traits.
                          nruns = 1,
                          big = FALSE,
                          drop.lod = 1.5, ## LOD drop to keep.
                          pheno.col = seq(n.phe), ## Traits used for hotspots.
                          ...)
{
  ## PHASE 1: Set up cross object. Needed in phases 2.
  ##
  ## Input files:
  ##       cross.RData: cross object and possibly lod.thrs or other objects.
  ##       params.txt
  ##
  ## Output files:
  ##       Phase1.RData
  ##

  ## Get any parameters in file. These will overwrite passed arguments.
  eval(parse(file.path(dirpath, params.file)))

  ## Cross object. Load if not done already.
  if(!exists(cross.name))
    load(file.path(dirpath, cross.file))
  
  ## Change name of cross object to "cross" for internal use.
  if(cross.name != "cross")
    cross <- get(cross.name)

  ## cross.index is used when multiple phase1 jobs are spawned.
  cross.index <- as.numeric(index)

  if(big) {
    ## Used for big data run.
    ## Each run is separate permutation.
    ## Make sure n.split is equal to n.perm. n.perm set to 1 in big.phase1.
    n.split <- max(1, n.perm)
    big.phase1(dirpath, cross.index, params.file, cross, lod.thrs, Nmax, n.perm,
               n.split, drop.lod = drop.lod, ...)
  }
  else {
    ## Used for studying properties of qtlhot. Not compatible with covariates.

    n.perm <- ceiling(n.perm / n.split)
      
    if(cross.index > 0 & nruns > 1) {
      ## Keep markers but re-simulate genotypes each time.
      mymap <- pull.map(cross)
      cross <- sim.cross(map = mymap, n.ind = n.ind, type = class(cross)[1])
      cross <- calc.genoprob(cross, step=0)
      
      ## Simulate new phenotypes.
      cross <- sim.null.pheno.data(cross, n.phe, latent.eff, res.var)
    }
    else
      cross.index <- 0

    Nmax <- min(Nmax, n.phe)
    
    ## Save all relevant objects for later phases.
    save(cross, n.phe, latent.eff, res.var, Nmax, n.perm, n.split,
         alpha.levels, lod.thrs, cross.index, big, pheno.col,
         file = file.path(dirpath, "Phase1.RData"),
         compress = TRUE)
  }

  ## Need to write a file with n.groups lines and group.size columns.
  ## The file groups.txt is examined by parallelizer (SOAR).
  write.table(c(cross.index, n.split),
              file = file.path(dirpath, "groups.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
####################################################################################
qtlhot.phase2 <- function(dirpath, index = NULL, ...,
                          ## Following are loaded with Phase1.RData created in big.phase1.
                          n.split, cross, Nmax, n.perm, lod.thrs, n.phe, alpha.levels,
                          drop.lod,
                          ##
                          batch.effect = NULL, addcovar = NULL, intcovar = NULL,
                          big = FALSE, verbose = FALSE)
{
  ## PHASE 2: NL,N and WW permutations. 
  ##          Slow. Run on condor nodes. Sized by n.perm.
  ##
  ## Steps for time calculations: n.perm * nind(cross) * nmar(cross) * nphe(cross)
  ## Input files:
  ##       Phase1.RData: cross, n.phe, n.perm, lod.thrs, 
  ##
  ## Output file (one per invocation):
  ##       permi.RData
  ##

  cross.index <- scan("groups.txt", 0)[1]
  ## Just to fool R build.
  pheno.col <- 0

  ## Load Phase 1 computations.
  infile <- "Phase1.RData"
  load(file.path(dirpath, infile))

  ## Quality check of index.
  if(missing(index))
    parallel.error(2, 2, index)
  index <- as.integer(index)
  if(is.na(index))
    parallel.error(3, 2, index)
  if(index < 1 | index > n.split)
    parallel.error(4, 2, index)

  if(big)
    return(big.phase2(dirpath, index, batch.effect = batch.effect,
                      addcovar = addcovar, intcovar = intcovar, ..., verbose = verbose))

  outfile <- paste("perm", ".", cross.index, "_", index, ".RData", sep = "")

  if(!is.null(batch.effect)) ## This may not be enough for some traits...
    covars <- sexbatch.covar(cross, batch.effect, verbose = TRUE)
  else
    covars <- list(addcovar = addcovar, intcovar = intcovar)
  
  ## Following is in hotperm stuff.
  ## filter.threshold is big loop. Have perm loop within it. n.phe*n.perm runs of scanone per dataset.
  ## For simulation study, have many datasets!

  ## Creates max.N of size n.perm x n.lod and max.lod.quant of size n.perm x Nmax.
  ## Size of n.perm determines the run time.
  mycat("hotperm", verbose)
  NLN <- hotperm(cross, Nmax, n.perm, alpha.levels, lod.thrs, drop.lod, verbose)

  mycat("scanone", verbose)
  scanmat <- scanone(cross, pheno.col = pheno.col, method = "hk",
                     addcovar=covars$addcovar, intcovar=covars$intcovar, ...)

  ## Reduce to high LOD scores.
  mycat("highlod", verbose)
  highobj <- highlod(scanmat, min(lod.thrs), drop.lod, restrict.lod = TRUE)
  rm(scanmat)
  gc()
  
  WW <- ww.perm(highobj, n.perm, alpha.levels, lod.thrs, n.perm)
  
  save(NLN, WW, index, lod.thrs, alpha.levels, drop.lod, covars,
       file = outfile, compress = TRUE)
}
####################################################################################
qtlhot.phase3 <- function(dirpath, index = NULL, ...,
                          ## Following are loaded with Phase1.RData created in big.phase1.
                          cross.index, n.split, n.perm, lod.thrs, Nmax, NLN, WW,
                          alpha.levels, cross, n.phe, latent.eff, res.var,
                          drop.lod,
                          ##
                          dirpath.save = dirpath, big = FALSE, verbose = FALSE)
{
  ## PHASE 3: Sample Markov chain (MCMC). Parallelize.
  ##          Fast: Run on scheduler.
  ##          Slow. Run on condor nodes.
  ##
  ## Input files:
  ##       Phase1.RData
  ##       permN.RData (multiple files with N = 1,...,n.split
  ##       cross.RData: cross object and possibly lod.thrs or other objects.
  ##
  ## Output files (one per invocation):
  ##       Phase3.RData: max.N, max.lod.quant
  ##
  ## See Phase 2 for explanation of NLNpermi.RData files.
  ## All NLNpermi.RData files are combined to make max.N, max.lod.quant.

  ## Just to fool R build.
  pheno.col <- 0
  covars <- NULL

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  if(big)
    return(big.phase3(dirpath, index, cross.index))

  ## This could be done once, but it would require splitting this phase in two.
  ## Besides, it is quite fast.
  ## Read in saved BIC scores and combine into one object.
  outfile <- paste("perm", ".", cross.index, "_", "*.RData", sep = "")
  filenames <- list.files(dirpath, outfile)
  if(!length(filenames))
    parallel.error(5, 3, index)
  if(length(filenames) != n.split)
    parallel.error(11, 3, index)

  n.perms <- n.perm * n.split
  n.lod <- length(lod.thrs)
  
  max.ww <- matrix(NA, n.perms, n.lod)
  class(max.ww) <- c("ww.perm", class(max.ww))
  qh.out <- list(max.lod.quant = matrix(NA, n.perms, Nmax),
                 max.N = max.ww)
  class(qh.out) <- c("hotperm", "list")

  i.perm <- seq(n.perm)
  for(i in seq(length(filenames))) {
    load(file.path(dirpath, filenames[i]))

    ## Do any quality checking here.
    
    qh.out$max.N[i.perm, ] <- NLN$max.N
    qh.out$max.lod.quant[i.perm, ] <- NLN$max.lod.quant
    max.ww[i.perm, ] <- WW
    
    i.perm <- i.perm + n.perm
  }
  attr(max.ww, "lod.thrs") <- lod.thrs
  attr(max.ww, "alpha.levels") <- alpha.levels
  attr(qh.out, "lod.thrs") <- lod.thrs
  attr(qh.out, "alpha.levels") <- alpha.levels

  qh.thrs <- summary(qh.out, alpha.levels)
  ww.thrs <- summary(max.ww, alpha.levels)

  ## Now compare permutations to original cross object from Phase1.

  out.sim <- filter.threshold(cross, pheno.col, latent.eff, res.var,
                              lod.thrs, drop.lod, seq(Nmax), n.perms, alpha.levels,
                              qh.thrs, ww.thrs,
                              addcovar = covars$addcovar, intcovar = covars$intcovar,
                              verbose = verbose, ...)

  outfile <- paste("Phase3", ifelse(cross.index > 0, cross.index, ""), ".RData", sep = "")
  
  save(latent.eff, Nmax, out.sim, 
       file = file.path(dirpath.save, outfile),
       compress = TRUE)
}
