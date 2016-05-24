#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# hard.thresh.scalar : adjust x according to a scalar threshold                #
#    Original code provided by Dr. Dan Yang, modified only in commentary.      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x   : matrix to be adjusted                                                #
#   thr : threshold value                                                      #
# Outputs                                                                      #
#    adjusted x                                                                #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
hard.thresh.scalar <- function(x, 
                               thr){
  ans <- x
  ans[x < thr & x > -thr]  <-  0
  ans
}
