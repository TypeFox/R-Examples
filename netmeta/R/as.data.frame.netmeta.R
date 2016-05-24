as.data.frame.netmeta <- function(x, row.names=NULL,
                                  optional=FALSE,
                                  details=FALSE, ...){
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")

  ## Remove element 'call' from object of class meta to get rid
  ## of an error message in meta-analyses with six studies:
  ## 'Error: evaluation nested too deeply: infinite recursion ...'
  ##
  ## NB: Element 'call' which is of length six contains information
  ##     on the function call.
  ##
  x$call <- NULL
  
  sel <- as.vector(lapply(x, length) == length(x$studlab))
  
  res <- as.data.frame(x[names(x)[sel]], ...)
  
  if (!details)
    res <- res[, !(names(res) %in% c("treat1.pos", "treat2.pos",
                                     "lower.nma.fixed", "upper.nma.fixed",
                                     "lower.nma.random", "upper.nma.random",
                                     "leverage.fixed"))]
  
  attr(res, "version") <- packageDescription("netmeta")$Version
  
  res
}
