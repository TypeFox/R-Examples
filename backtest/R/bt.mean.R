################################################################################
##
## $Id: bt.mean.R 1300 2008-08-27 21:01:11Z zhao $
##
## Calculates the means by column and returns a
## 1 x length(x) array where the values are the means
## of the columns.
##
################################################################################

.bt.mean <- function(x){
  
  stopifnot(
            is.array(x)
            )

  x <- array(colMeans(x, na.rm = TRUE), dim = c(1, ncol(x)),
             dimnames = list("MEAN", dimnames(x)[[2]]))

  x
}

