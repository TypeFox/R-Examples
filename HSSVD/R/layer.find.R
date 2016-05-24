#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# layer.find Call to FIT-SSVD methods developed in arXiv:1112.2433             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# By default, FIT-SSVD uses hard thresholding, values for which are obtained   #
#   using the median bootstrap methodology (Algorithm 3)                       #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x       : matrix of observed data                                          #
#   u0      : matrix of left singular vectors                                  #
#   v0      : matrix of right singular vectors                                 #
#  arg.list : inputs for FIT.SSVD                                              #
# Outputs                                                                      #
#    returns object created by SSVD.iter.thresh, which includes                #
#  u         : estimator for left singular vectors                             #
#  v         : estimator for right singular vectors                            #
#  d         : estimator for singular values                                   #
#  k.u       : number of non-zero elements in each vector of u                 #
#  k.v       : number of non-zero elements in each vector of v                 #
#  n.step    : number of iterations required                                   #
#  sigma.hat : estimated standard deviation                                    #
#  thr.u     : final threshold for u                                           #
#  thr.v     : final threshold for v                                           #
#  dist      : maximum distance between subspaces                              #
#  messages  : information return by internal routines                         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

layer.find <- function(x, 
                       u0, 
                       v0,
                       arg.list){

  arg.list[['x']] <- x
  arg.list[['u.old']] <- u0
  arg.list[['v.old']] <- v0

  result <- do.call(SSVD.iter.thresh,arg.list)

  return(result)
}
