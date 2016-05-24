##
## This file contains ORIGINAL COPIES of methods from the base package stats (2014-10-31)
## and the extension packages car (2014-10-10) and multcomp (2013-10-04).
##



## Function copied from package stats to make anova.lm work.
## using qr(<lmm>)  as interface to  <lmm>$qr :
qr.lmm <- function(x, ...) {
  if(is.null(r <- x$qr))
    stop("lm object does not have a proper 'qr' component.
 Rank zero or should not have used lm(.., qr=FALSE).")
  r
}
