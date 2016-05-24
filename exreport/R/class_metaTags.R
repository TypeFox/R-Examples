is.metaTags <- function(x) {
  is(x, "metaTags")
}

.metaTags <- function(
  alias="",         # A name for a given operation, used in rendering.
  alpha="",         # Significance Level of a statistical test.
  context="",       # Name of the source experiment for procedure.
  distribution="",  # Name of the  distribution of a statistical test.
  method="" ,       # (subjet to change) A parametrized method of a procedure. 
  objetive="",      # (subjet to change) objetive function of a procedure.
  outcome="",       # Interpretation of the testing procedure.
  pvalue="",        # Statistical value for several tests.
  scope="",         # (subjet to change) A parametrized method of a procedure.
  statistic="",     # Statistical value computed for a test.
  target="",        # Target variable of a procedure.
  title=""          # Name of a reportable object.
  ) {
  
  fargs <- as.list(environment())
  
  error <- lapply(fargs, function(x) !is.null(x))
  fargs <- fargs[unlist(error)]
  
  cond <- lapply(fargs, function(x) x!="")
  fargs <- fargs[unlist(cond)]
  
  class(fargs) <- c("metaTags")
  fargs
}


.updateTags <- function(t1, t2) {
  for (name in names(t2)) t1[[name]] = t2[[name]]
  t1
}