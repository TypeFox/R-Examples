###############################################################
#This function is a wrapper for the apply function            #
#that has a progress bar - Author: Mark Heckman               #          #                                                             #
#http://ryouready.wordpress.com/2010/01/11/                   #
#       progress-bars-in-r-#part-ii-a-wrapper-for-apply-      #
#       functions/                                            #
###############################################################

apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,
                       style = 3)

  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                           curVal +1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}