###########################################################################/**
# @set "class=environment"
# @RdocMethod objectSize
#
# @title "Gets the size of an environment in bytes"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{envir}{An @see "base::environment".}
#   \item{...}{Arguments passed to @see "base::ls".}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   Internally @see "utils::object.size" is used.
# }
#
# \keyword{attribute}
# \keyword{utilities}
#*/###########################################################################
setMethodS3("objectSize", "environment", function(envir, ...) {
  ## Keep track of already scanned environments
  ## in order to avoid endless recursion.
  args <- list(...)
  .scannedEnvs <- args$.scannedEnvs
  if (is.null(.scannedEnvs)) .scannedEnvs <- list()

  alreadyScanned <- function(envir) {
    if (!is.environment(envir)) return(FALSE)
    for (env in .scannedEnvs) {
      if (identical(env, envir)) return(TRUE)
    }
    FALSE
  }

  ## Get all objects in the environment
  args <- list(envir=envir, all.names=TRUE, ...)
  args$.scannedEnvs <- NULL
  names <- do.call(ls, args=args)
  # Nothing to do?
  if (length(names) == 0L) return(0)

  ## Avoid scanning the current environment again
  .scannedEnvs <- c(.scannedEnvs, envir)

  size <- 0
  for (name in names) {
    obj <- envir[[name]]
    if (!alreadyScanned(obj)) {
      size <- size + objectSize(obj, .scannedEnvs=.scannedEnvs)
    }
  }

  size
})




############################################################################
# HISTORY:
# 2015-01-27
# o BUG FIX: objectSize() for environment could result in infinite
#   recursive calls if there circular dependencies between environments.
# 2009-10-26
# o Added objectSize() for environments.
############################################################################
