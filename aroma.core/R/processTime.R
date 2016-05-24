###########################################################################/**
# @RdocDefault processTime
#
# @title "Gets the running time of the R process and its children processes"
#
# \description{
#   @get "title".
#   This function is a safe wrapper for @see "base::proc.time", which might
#   not exist on all platforms.  It "determines how much time (in seconds) 
#   the currently running \R process already consumed".  In addition it adds
#   descriptive names of the returned values.  
#   For more details, see @see "base::proc.time".
# }
#
# @synopsis
#
# \arguments{
#  \item{since}{An optional @numeric @vector to be subtracted from the value
#    of @see "base::proc.time".  This is useful for calculating "lap times".}
#  \item{units}{A @character string specifying the unit of the 
#    returned values.}
#  \item{fmtstr}{If given, a format string to convert the times to strings
#    via @see "base::sprintf".}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @numeric @vector of length 5.
#   For more details, see @see "base::proc.time".
# }
#
# @author
#
# \seealso{
#   @see "base::proc.time".
#   @see "base::system.time".
#   @see "base::gc.time".
# }
#
# @keyword "programming"
# @visibility "private"
#*/###########################################################################
setMethodS3("processTime", "default", function(since=NULL, units=c("seconds", "milliseconds", "minutes", "hours", "days"), fmtstr=NULL, ...) {
  if (!exists("proc.time", mode="function"))
    return(rep(as.double(NA), times=5L));

  units <- match.arg(units);

  time <- proc.time();

  if (!is.null(since))
    time <- time - since;

  if (units != "seconds") {
    if (units == "milliseconds") {
      time <- 1000 * time;
    } else if (units == "minutes") {
      time <- time / 60;
    } else if (units == "hours") {
      time <- time / 3600;
    } else if (units == "days") {
      time <- time / (86400);
    }
  }

  if (!is.null(fmtstr)) {
    time <- sprintf(fmtstr, time);
  }

  names(time) <- c("user", "system", "total", "userChilds", "systemChilds");

  time;
}, private=TRUE); # processTime()


############################################################################
# HISTORY:
# 2006-03-09
# o Created.
############################################################################
