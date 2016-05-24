# CHNOSZ/util.program.R
# various programming-related functions

caller.name <- function(n=2) {
  # returns the name of the calling function n frames up
  # (n=2: the caller of the function that calls this one)
  # or character() if called interactively
  if(sys.nframe() < n) name <- character()
  else {
    sc <- sys.call(-n)[[1]]
    name <- try(as.character(sc),silent=TRUE)
    # also return character() if the value from sys.call is
    # the function itself (why does this sometimes happen,
    # e.g. when called from affinity()?)
    if(class(name)=="try-error") name <- character()
  }
  return(name)
}

palply <- function(varlist, X, FUN, ...) {
  # a wrapper function to run parLapply if length(X) >= thermo$opt$paramin
  # and package 'parallel' is available, otherwise run lapply
  if(length(X) >= get("thermo")$opt$paramin) {
    # Use option mc.cores to choose an appropriate cluster size.
    # and set max at 2 for now (per CRAN policies)
    nCores <- min(getOption("mc.cores"), 2)
    # don't load methods package, to save startup time - ?makeCluster
    cl <- parallel::makeCluster(nCores, methods=FALSE)
    # export the variables and notify the user
    if(! "" %in% varlist) {
      parallel::clusterExport(cl, varlist)
      msgout(paste("palply:", caller.name(4), "running", length(X), "calculations on",
        nCores, "cores with variable(s)", paste(varlist, collapse=", "), "\n"))
    } else {
      msgout(paste("palply:", caller.name(4), "running", length(X), "calculations on",
        nCores, "cores\n"))
    }
    # run the calculations
    out <- parallel::parLapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
  } else out <- lapply(X, FUN, ...)
  return(out)
}
