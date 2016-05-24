#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# soft.thresh.scalar : adjust x according to a scalar threshold                #
#    Original code provided by Dr. Dan Yang, modified in commentary added      #
#    explicit selection criteria.                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x   : matrix to be adjusted                                                #
#   thr : threshold value                                                      #
# Outputs                                                                      #
#    adjusted x                                                                #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
soft.thresh.scalar <- function(x, 
                               thr){
  ans <- x
  ans[x > thr]  <-  x[x > thr] - thr
  ans[x < -thr]  <-  x[x < -thr] + thr
  ans[x < thr & x > -thr]  <-  0
  ans
}
