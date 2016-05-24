#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# soft.thresh : initiates hard thresholding                                    #
#    Original code provided by Dr. Dan Yang, modified only in commentary.      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x   : matrix to be adjusted                                                #
#   thr : threshold value                                                      #
# Outputs                                                                      #
#    adjusted x                                                                #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
soft.thresh <- function(x, 
                        thr){

  if(length(thr)==1) {
    #----------------------------------------------------------------------#
    # if single threshold value, apply to all elements of x                #
    #----------------------------------------------------------------------#
    return(soft.thresh.scalar(x,thr))
  } else {
    #----------------------------------------------------------------------#
    # if multiple threshold values, ensure dim matches matrix              #
    #----------------------------------------------------------------------#
    if(ncol(x) != length(thr)) {
      stop("soft.thresh: dimensions do not match")
    }
    #----------------------------------------------------------------------#
    # if multiple threshold values, apply each threshold to its column in x#
    #----------------------------------------------------------------------#
    ans <- x
    for(i in 1:length(thr)){
      ans[,i] <- soft.thresh.scalar(x[,i],thr[i])
    }
    return(ans)
  }
}
