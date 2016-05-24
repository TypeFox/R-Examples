###########################################################################/**
# @RdocClass RspShSourceCode
#
# @title "The RspShSourceCode class"
#
# \description{
#  @classhierarchy
#
#  An RspShSourceCode object is an @see "RspSourceCode" holding R source code.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{@character strings.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspShSourceCode", function(...) {
  extend(RspSourceCode(...), "RspShSourceCode");
})



#########################################################################/**
# @RdocMethod evaluate
# @aliasmethod findProcessor
#
# @title "Evaluates the shell (sh) code"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{envir}{The @environment in which the RSP string is evaluated.}
#   \item{args}{A named @list of arguments assigned to the environment
#     in which the RSP string is parsed and evaluated.
#     See @see "R.utils::cmdArgs".}
#   \item{output}{A @character string specifying how the RSP output
#     should be handled/returned.}
#   \item{...}{Optional arguments passed to @see "base::eval".}
# }
#
# \value{
#  If \code{output="stdout"}, then @NULL is returned and the RSP output
#  is sent to the standard output.
#  Note that this is output is "buffered", meaning it will be sent to
#  standard output upon completion.  This is a limitation of R.
#  If \code{output="RspStringProduct"}, then the output is captured
#  and returned as an @see "RspStringProduct" with attributes set.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("evaluate", "RspShSourceCode", function(object, envir=parent.frame(), args="*", output=c("RspStringProduct", "stdout"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  evalSh <- function(text, ...) {
    pathnameT <- tempfile(pattern="RSP-sh-", fileext=".sh")
    writeLines(text, con=pathnameT);
    on.exit({
      file.remove(pathnameT);
    })
    res <- system2("sh", args=list(pathnameT), stdout=TRUE);
    res <- paste(res, collapse="\n");
    res;
  } # evalSh()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'args':
  args <- cmdArgs(args=args);

  # Argument 'output':
  output <- match.arg(output);


  code <- object;

  # Assign arguments to the parse/evaluation environment
  attachLocally(args, envir=envir);

  # Evaluate R source code and capture output
  res <- evalSh(code);

  if (output == "RspStringProduct") {
    res <- RspStringProduct(res, type=getType(object));
  } else {
    cat(res);
    res <- NULL;
  }

  res;
}, createGeneric=FALSE) # evaluate()


setMethodS3("findProcessor", "RspShSourceCode", function(object, ...) {
  function(...) {
    evaluate(...);
  }
}) # findProcess()



##############################################################################
# HISTORY:
# 2013-08-04
# o Added argument 'output' to evaluate() for RspShSourceCode.
# 2013-03-14
# o Created from RspRSourceCode.
##############################################################################
