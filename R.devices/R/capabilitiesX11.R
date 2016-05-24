###########################################################################/**
# @RdocFunction capabilitiesX11
#
# @title "Check whether current R session supports X11 or not"
#
# \description{
#   @get "title".
#
#   Contrary to \code{capabilities("X11")} which only checks for X11
#   support on startup [1], this function checks whether X11 is supported
#   when it is called.  This is done by querying a temporary R session.
# }
#
# @synopsis
#
# \arguments{
#  \item{reset}{If @TRUE, any previously obtained results are ignored,
#   otherwise not.}
#  \item{...}{(optional) @character strings of command-line options
#   to \code{Rscript}.}
# }
#
# \value{
#   Returns @TRUE if X11 is supported, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#  @see base::capabilities
# }
#
# \references{
#  [1] R-devel thread 'capabilities("X11"): Force refresh from within R?
#      (... and minor documentation issue)' on 2015-05-06,
#      \url{https://stat.ethz.ch/pipermail/r-devel/2015-May/071100.html}\cr
# }
#
# @keyword device
# @keyword internal
#*/###########################################################################
capabilitiesX11 <- local({
  capableOfX11 <- NA
  function(reset=FALSE, ...) {
    if (is.na(capableOfX11) || reset) {
      bin <- file.path(R.home("bin"), "Rscript")
      cmd <- "cat(capabilities('X11'))"
      args <- c(..., "-e", dQuote(cmd))
      value <- system2(bin, args=args, stdout=TRUE)
      capableOfX11 <<- isTRUE(as.logical(value))
    }
    capableOfX11
  }
})
