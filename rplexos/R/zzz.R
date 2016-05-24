.onLoad <- function(libname, pkgname) {
  # Set options
  op <- options()
  op.rplexos <- list(
    rplexos.tiebreak = "last",
    rplexos.debug = FALSE
  )
  toset <- !(names(op.rplexos) %in% names(op))
  if(any(toset)) options(op.rplexos[toset])
  
  # By default, turn off parallel queries
  assign("rplexos_globals", new.env(), envir=parent.env(environment()))
  if (!"cluster" %in% ls(rplexos_globals))
    assign("cluster", NULL, rplexos_globals)
  start_parallel_rplexos(1, TRUE)

  invisible()
}
