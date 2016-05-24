## Directory path is indicated with "dirpath" character string
##
## I/O files are in internal R format (RData) using load/save commands
## except "groups.txt", which is used to set width of parallelization.
##
## I/O files are on AFS in
## /u/y/a/yandell/public/html/sysgen/qtlnet/condor
## or equivalently via URL
## http://www.stat.wisc.edu/~yandell/sysgen/qtlnet/condor
##
####################################################################################
parallel.qtlnet <- function(phase, index = 1, ..., dirpath = ".")
{
  switch(phase,
         qtlnet.phase1(dirpath, index, ...),
         qtlnet.phase2(dirpath, index, ...),
         qtlnet.phase3(dirpath, index, ...),
         qtlnet.phase4(dirpath, index, ...),
         qtlnet.phase5(dirpath, index, ...),
         parallel.error(1, phase, index))
  
  parallel.error(0, phase, index)
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
                     "phase number not between 1 and 5",
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
qtlnet.phase1 <- function(dirpath, index = NULL,
                          params.file = "params.txt",
                          cross.file = "cross.RData",
                          cross.name = "cross",
                          pheno.col = 1:13,
                          max.parents = 12,
                          threshold = 3.83,
                          nruns = 100,
                          step = 1,
                          addcov = NULL, intcov = NULL,
                          nSamples = 1000, thinning = 20,
                          n.groups = NULL, group.size = 50000,
                          M0 = NULL, init.edges = NULL,
                          ...)
{
  ## PHASE 1: Initiation. Needed in phases 2 and 3.
  ##          Fast: Run on scheduler.
  ##
  ## Input files:
  ##       cross.RData
  ##       params.txt
  ##
  ## Output files:
  ##       Phase1.RData
  ##       groups.txt
  ##
  ## The "groups.txt" file is used to determine Phase2 width.

  ## Get any parameters in file. These will overwrite passed arguments.
  eval(parse(file.path(dirpath, params.file)))

  ## Cross object. Load if not done already.
  if(!exists(cross.name))
    load(file.path(dirpath, cross.file))

  ## Change name of cross object to "cross" for internal use.
  if(cross.name != "cross")
    cross <- get(cross.name)

  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]])))
    cross <- calc.genoprob(cross, step = step)
  
  ## Break up into groups to run on several machines.
  ## 54 groups of ~1000, for a total of 53248 scanone runs.
  parents <- parents.qtlnet(pheno.col, max.parents)
  groups <- group.qtlnet(parents = parents, n.groups = n.groups,
                         group.size = group.size)
  
  ## Save all relevant objects for later phases.
  save(cross, pheno.col, max.parents, threshold, parents, groups, nruns,
       addcov, intcov, nSamples, thinning, M0, init.edges,
       file = file.path(dirpath, "Phase1.RData"),
       compress = TRUE)
  
  ## Need to write a file with n.groups lines and group.size columns.
  write.table(groups,
              file = file.path(dirpath, "groups.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
####################################################################################
qtlnet.phase2 <- function(dirpath, index = NULL, ...,
                          groups, cross, pheno.col, threshold, max.parents, parents,
                          verbose = FALSE)
{
  ## PHASE 2: Compute BIC scores. Parallelize.
  ##          Slow. Run on condor nodes.
  ##
  ## Input files:
  ##       Phase1.RData
  ##
  ## Output file (one per invocation):
  ##       bicN.RData (with N = index)
  ##
  ## Argument:
  ##       index (1 to number of lines in groups.txt file)
  ##
  ## The "groups.txt" file (created in Phase1) is used to determine Phase 2 width.
  ## That is, groups.txt has 54 lines, hence 54 separate condor runs
  ## to produce bic1.RData through bic54.RData.
  
  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  ## Quality check of index.
  if(missing(index))
    parallel.error(2, 2, index)
  index <- as.integer(index)
  if(is.na(index))
    parallel.error(3, 2, index)
  if(index < 1 | index > nrow(groups))
    parallel.error(4, 2, index)

  ## Pre-compute BIC scores for selected parents.
  bic <- bic.qtlnet(cross, pheno.col, threshold,
                    max.parents = max.parents,
                    parents = parents[seq(groups[index,1], groups[index,2])],
                    verbose = verbose,
                    ...)
  
  save(bic,
       file = file.path(dirpath, paste("bic", index, ".RData", sep = "")),
       compress = TRUE)
}
####################################################################################
qtlnet.phase3 <- function(dirpath, index = NULL, ...,
                          groups, bic, cross, pheno.col, max.parents)
{
  ## PHASE 3: Sample Markov chain (MCMC). Parallelize.
  ##          Fast: Run on scheduler.
  ##          Slow. Run on condor nodes.
  ##
  ## Input files:
  ##       Phase1.RData
  ##       bicN.RData (multiple files with N = 1,...,nrow(groups))
  ##
  ## Output files (one per invocation):
  ##       Phase3.RData
  ##
  ## See Phase 2 for explanation of bicN.RData files.
  ## All bicN.RData files are combined to make saved.scores.

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  ## This could be done once, but it would require splitting this phase in two.
  ## Besides, it is quite fast.
  ## Read in saved BIC scores and combine into one object.
  filenames <- list.files(dirpath, "bic.*RData")
  if(!length(filenames))
    parallel.error(5, 3, index)
  if(length(filenames) != nrow(groups))
    parallel.error(11, 3, index)
  
  bic.group <- list()
  for(i in seq(length(filenames))) {
    load(file.path(dirpath, filenames[i]))

    ## Do any quality checking here.
    
    bic.group[[i]] <- bic
  }
  
  saved.scores <- bic.join(cross, pheno.col, bic.group, max.parents = max.parents)
  save(saved.scores,
       file = file.path(dirpath, "Phase3.RData"),
       compress = TRUE)
}
####################################################################################
qtlnet.phase4 <- function(dirpath, index = NULL, ...,
                          nruns, cross, pheno.col, threshold, max.parents, saved.scores,
                          init.edges, M0, addcov, intcov, nSamples, thinning)
{
  ## PHASE 4: Sample Markov chain (MCMC). Parallelize.
  ##          Slow. Run on condor nodes.
  ##
  ## Input files:
  ##       Phase1.RData
  ##       Phase3.RData
  ##
  ## Output files (one per invocation):
  ##       mcmcN.RData (actual names will use temporary strings for XXX)

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  ## Load Phase 3 computations (saved.scores).
  load(file.path(dirpath, "Phase3.RData"))

  ## Quality check of index.
  if(missing(index))
    parallel.error(6, 4, index)
  index <- as.integer(index)
  if(is.na(index))
    parallel.error(7, 4, index)
  if(index < 1 | index > nruns)
    parallel.error(8, 4, index)

  ## Run MCMC with randomized initial network.
  mcmc <- mcmc.qtlnet(cross, pheno.col, threshold = threshold,
                      max.parents = max.parents,
                      saved.scores = saved.scores,
                      init.edges = init.edges, M0 = M0,
                      addcov = addcov, intcov = intcov,
                      nSamples = nSamples, thinning = thinning,
                      ...)
  save(mcmc,
       file = file.path(dirpath, paste("mcmc", index, ".RData", sep = "")),
       compress = TRUE)
}
####################################################################################
qtlnet.phase5 <- function(dirpath, index = NULL, missing.ok = FALSE, ...,
                          nruns,
                          verbose = FALSE)
{
  ## PHASE 5: Combine results for post-processing.
  ##          Fast: Run on scheduler.
  ##
  ## Input files (results of Phase 3):
  ##       Phase1.RData
  ##       mcmcN.RData (N = 1, ..., nruns)
  ##
  ## Output files:
  ##       Final.RData
  ##
  ## All "mcmc.*RData" files in "dirpath" are combined to creat output file.

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  ## Combine outputs together.
  filenames <- list.files(dirpath, "mcmc.*RData")
  if(!length(filenames))
    parallel.error(9, 5, index)
  if(length(filenames) != nruns & !missing.ok)
    parallel.error(10, 5, index)
  
  for(i in seq(length(filenames))) {
    if(verbose)
      print(i)
    load(file.path(dirpath, filenames[i]))
    mcmc <- legacy.qtlnet(mcmc, codes = (i == 1))
    if(i == 1)
      out.qtlnet <- mcmc
    else
      out.qtlnet <- c(out.qtlnet, mcmc)
  }

  save(out.qtlnet,
       file = file.path(dirpath, "Final.RData"),
       compress = TRUE)
}
