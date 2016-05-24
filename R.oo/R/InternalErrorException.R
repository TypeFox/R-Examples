###########################################################################/**
# @RdocClass InternalErrorException
#
# @title "InternalErrorException represents internal errors"
#
# \description{
#  @classhierarchy
#  
#  @get "title" that are likely to be due to implementation errors done by
#  the author of a specific package and not because the user made an error.
#  Errors that are due to unexpected input to functions etc falls under
#  this error type.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Any arguments accepted by @see "Exception"}.
#   \item{package}{The name (@character string) of the package where the 
#     error exists. Can also be a @see "Package" object. If @NULL, the
#     source of the error is assumed to be unknown.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# \seealso{
#   For detailed information about exceptions see @see "Exception".
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setConstructorS3("InternalErrorException", function(..., package=NULL) {
  if (!is.null(package) && !inherits(package, "Package")) {
    package <- Package(as.character(package));
  }

  extend(Exception(...), "InternalErrorException",
    .package=package
  )
})


###########################################################################/**
# @RdocMethod getPackage
#
# @title "Gets the suspicious package likely to contain an error"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "Package" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getPackage", "InternalErrorException", function(this, ...) {
  this$.package;
})



###########################################################################/**
# @RdocMethod getMessage
#
# @title "Gets the message of the exception"
#
# \description{
#  @get "title" and adds a message that clarifies that the error is likely
#  due to an internal error and not due to the user. It also gives information
#  how to contact the maintainer or author of the suspicous package. This 
#  information is retrieved from the DESCRIPTION file of the package. To help
#  the package developer, information about the current version of R, the
#  current version of the package etc are also returned so the user can 
#  easily cut'n'paste that information into a bug report.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
# \keyword{error}
#*/###########################################################################
setMethodS3("getMessage", "InternalErrorException", function(this, ...) {
  msg <- getMessage.Exception(this);
  msg <- paste(msg, " This error is likely to be due to an internal error", sep="");

  pkg <- getPackage(this);
  if (!is.null(pkg)) {
    msg <- paste(msg, " related to package ", getName(pkg), " v", getVersion(pkg), ". Please report this problem to the maintainer ", getMaintainer(pkg), " or the author ", getAuthor(pkg), " of that package", sep="");
  }

  R.oo <- Package("R.oo");
  msg <- paste(msg, ". Do not forget to report that you are using R v", getVersion(Package("base")), " on a ", R.Version()$platform, " platform together with R.oo v", getVersion(R.oo), ".", sep="");
  msg;
})


############################################################################
# HISTORY:
# 2012-11-13
# o ROBUSTNESS: Now getMessage() for InternalErrorException setups an
#   Package("R.oo") object, instead of assuming that R.oo is loaded.
# 2007-04-07
# o Removed reportBug() for InternalErrorException.  Never completed/used.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2003-12-02
# o Bug report now generates a form process by the braju.com server.
# 2003-04-15
# o First trial version of a bug report system from within R that 
#   automatically fills in the version information since that is the most
#   commonly forgotten information.
# o Created.
############################################################################
