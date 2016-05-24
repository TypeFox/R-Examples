exportDoMPI <- function(cl, varlist, envir=.GlobalEnv) {
  genvir <- new.env(parent=emptyenv())

  if (missing(varlist)) {
    # Default to all functions from the specified environment
    for (nm in ls(envir, all.names=FALSE)) {
      obj <- get(nm, pos=envir, inherits=TRUE)
      if (is.function(obj)) {
        assign(nm, obj, pos=genvir)
      }
    }
  } else {
    for (nm in varlist) {
      assign(nm, get(nm, pos=envir, inherits=TRUE), pos=genvir)
    }
  }

  # Piggy backing only right now
  for (i in seq(length=clusterSize(cl))) {
    sendToWorker(cl, i,
       list(job=genvir, joblen=0, globaljob=TRUE, jobcomplete=TRUE,
            numtasks=0))
  }

  invisible(NULL)
}
