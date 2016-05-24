logistf.control <-
function(maxit=25, maxhs=5, maxstep=5, lconv=0.00001, gconv=0.00001, xconv=0.00001, collapse=TRUE){
  list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, gconv=gconv, xconv=xconv, collapse=collapse)
}

