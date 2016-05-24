## covr: skip=all

###########################################################################/**
# @RdocDefault abort
# @alias abort.condition
#
# @title "Aborts the current expression call"
#
# \description{
#  @get "title" and returns to the top level prompt/browser
#  \emph{without signalling a condition}.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{(optional) Objects coerced to @character and pasted together without a separator, or a @condition object. If no object are given, no message is printed.}
#   \item{call.}{If @TRUE, the call is added to the message, otherwise not.}
#   \item{domain}{Used to translate the message (see @see "base::gettext"). If @NA, messages will not be translated.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \details{
#   There are still cases where one can "catch" the abort and undo it, cf. [1].
# }
#
# \examples{\dontrun{
#  @include "../incl/abort.Rex"
# }}
#
# @author
#
# \seealso{
#   @see "throw".
#   @see "base::stop".
#   Internally, @see "base::invokeRestart"\code{("abort")} is utilized.
# }
#
# \references{
#   [1] R-devel thread '', Sept 11, 2012,
#   \url{https://stat.ethz.ch/pipermail/r-devel/2012-September/064838.html}.\cr
# }
#
# @keyword error
# @keyword internal
#*/###########################################################################
setMethodS3("abort", "condition", function(cond, ..., call.=TRUE, domain=NULL) {
  message <- conditionMessage(cond);
  call <- conditionCall(cond);
  if (is.null(call)) {
    msg <- sprintf("%s", .makeMessage("Abort", domain=domain));
  } else {
    call <- deparse(call);
    msg <- sprintf("%s %s", .makeMessage("Abort in", domain=domain), call);
  }
  msg <- sprintf("%s: %s\n", msg, message);
  cat(msg, file=stderr());
  abort();
})

setMethodS3("abort", "default", function(..., call.=TRUE, domain=NULL) {
  args <- list(...);
  if (nargs() > 0) {
    message <- .makeMessage(..., domain=domain);
    nframe <- sys.nframe();
    if (nframe <= 2) call. <- FALSE;
    if (call.) {
      call <- sys.call(which=nframe-2L);
      if (is.null(call)) {
        msg <- sprintf("%s", .makeMessage("Abort", domain=domain));
      } else {
        call <- deparse(call);
        msg <- sprintf("%s %s", .makeMessage("Abort in", domain=domain), call);
      }
      msg <- sprintf("%s: %s\n", msg, message);
    } else {
      msg <- sprintf("%s: %s\n", .makeMessage("Abort", domain=domain), message);
    }
    cat(msg, file=stderr());
  }

  # Now abort R.
  invokeRestart("abort");
})



############################################################################
# HISTORY:
# 2012-09-11
# o Now abort() immitates how stop() works but without the signalling
#   of a condition.
# 2012-09-10
# o ROBUSTNESS/CRAN POLICY: Updated abort() for condition to utilize
#   invokeRestart("abort").  This avoids having to call
#   .Internal(.signalCondition(...)).  It also means that the message
#   outputted by abort() no longer starts with a "Error in ...:" line.
# 2012-03-05
# o The abort() method is hidden and is not used by any R.oo methods.
#   Will keep it until it is fully certain that throw() for Exception
#   will work as expected without it.
# 2012-02-29
# o KNOWN ISSUES: abort() for 'condition' still uses .Internal().
# o Added abort(), which is available as default function as well as
#   a method for 'condition' objects.
# o Created.
############################################################################
