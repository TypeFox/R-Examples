###########################################################################/**
# @class "environment"
# @RdocMethod getName
#
# @title "Gets the name of an environment"
#
# \description{
#  @get "title", e.g. \code{"R_GlobalEnv"} or \code{"0x01ddd060"}.
# }
#
# @synopsis
#
# \arguments{
#   \item{env}{An @environment.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# \examples{
#   name <- getName(globalenv())
#   print(name)
#   stopifnot(identical(name, "R_GlobalEnv"))
#
#   getName(new.env())
# }
#
# @author
#
# \seealso{
#   \code{\link[base:environment]{environmentName}()}.
# }
#
# \keyword{programming}
#*/###########################################################################
setMethodS3("getName", "environment", function(env, ...) {
  # base::environmentName() was added to R v2.5.0
  if (exists("environmentName", mode="function")) {
    name <- environmentName(env);
  } else {
    name <- "";
  }

  if (name == "") {
    name <- capture.output(print.default(env));
    name <- name[1]; # Just in case
    name <- gsub("[<]*environment:[ ]*([^>]*)[>]", "\\1", name);
  }

  name;
})



############################################################################
# HISTORY:
# 2008-03-25
# o Added getName() for 'environment':s. It extends base::environmentName()
#   to return the "pointer" if not a name.  It is used by 
#   getInternalAddress() of Object in R.oo.
# o Created.
############################################################################
