###########################################################################/**
# @RdocFunction withPar
#
# @title "Evaluate an R expression with graphical parameters set temporarily"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{The R expression to be evaluated.}
#   \item{...}{Named options to be used.}
#   \item{args}{(optional) Additional named options specified as a named @list.}
#   \item{envir}{The @environment in which the expression should be evaluated.}
# }
#
# \value{
#  Returns the results of the expression evaluated.
# }
#
# \details{
#   Upon exit (also on errors), this function will reset \emph{all}
#   (modifiable) graphical parameters to the state of options available
#   upon entry.  This means any parameters \emph{modified} from evaluating
#   \code{expr} will also be undone upon exit.
# }
#
# @author
#
# @examples "../incl/withPar.Rex"
#
# \seealso{
#   Internally, @see "base::eval" is used to evaluate the expression,
#   and @see "graphics::par" to set graphical parameters.
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
withPar <- function(expr, ..., args=list(), envir=parent.frame()) {
  # Argument '.expr':
  expr <- substitute(expr)

  # Argument 'args':
  if (!is.list(args)) {
    throw("Argument 'args' is not a list: ", class(args)[1L])
  }

  # Argument 'envir':
  if (!is.environment(envir)) {
    throw("Argument 'envir' is not a list: ", class(envir)[1L])
  }

  # All arguments specified
  new <- c(list(...), args)

  # Set parameters temporarily (restore *all* upon exit)
  prev <- par(no.readonly=TRUE)
  on.exit(par(prev))
  if (length(new) > 0L) par(new)

  eval(expr, envir=envir)
} # withPar()


############################################################################
# HISTORY:
# 2014-09-15
# o ROBUSTNESS: Now withPar() also resets the the graphical parameters
#   at entry even if no explicit ones were specified.  This covers the
#   case when the 'expr' changes the parameters, e.g.
#   withPar({ par(cex=2L); plot(1:10) }).
# 2014-05-01
# o Created.
############################################################################
