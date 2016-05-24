#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# error.est.full : threshold level function: Algorithm 3 of FIT-SSVD described #
#    in  arXiv:1112.2433. Original code provided by Dr. Dan Yang, modified     #
#    only in commentary and treatment of messages.                             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x            : p x q matrix of observed data                               #
#   u            : p x k matrix of left singular vectors previously estimated  #
#   v            : q x k matrix of right singular vectors previously estimated #
#   n.err        : number of boostraps to be used (Default=100)                #
#   error.median : use median/mean                                             #
# Outputs                                                                      #
#   gaga         : threshold values                                            #
#   messages     : informative messages                                        #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
error.est.full <- function(x, 
                           u, 
                           v, 
                           n.err=100, 
                           error.median=FALSE){

  messages <- list()

  x <- as.matrix(x)
  u <- as.matrix(u)
  v <- as.matrix(v)

  p <- nrow(x)
  q <- ncol(x)
  r <- ncol(u)

  #--------------------------------------------------------------------------#
  # index set of rows which have all zero entries                            #
  #--------------------------------------------------------------------------#
  sel.u <- apply(u==0, 1, all)
  sel.v <- apply(v==0, 1, all)

  #--------------------------------------------------------------------------#
  # number of all zero rows in u & v                                         #
  #--------------------------------------------------------------------------#
  n.sel.u <- sum(sel.u)
  n.sel.v <- sum(sel.v)

  #--------------------------------------------------------------------------#
  # Size of the selected matrix                                              #
  #--------------------------------------------------------------------------#
  n.sel <- n.sel.u * n.sel.v

  if(n.sel < (log(p*(q-n.sel.v)) * p*(q-n.sel.v))){
    #----------------------------------------------------------------------#
    # If size is 'small', use central limit theorem                        #
    #----------------------------------------------------------------------#
    messages <- c(messages, "too many nonzeros central limit theorem used")
    gaga <- rep(mad(as.vector(x)) * sqrt(2*log(p)), r)
  } else{
    #----------------------------------------------------------------------#
    # If size is large, use median/mean of M bootstrap replications        #
    #----------------------------------------------------------------------#
    samp <- sample(x=n.sel, size=(q-n.sel.v) * p * n.err, replace = TRUE)
    ans <- abs(matrix(x[sel.u, sel.v][samp], ncol=q-n.sel.v) %*% v[!sel.v,])
    if(error.median==FALSE){
      gaga <- apply(ans, 2, function(x){mean(apply(matrix(x,n.err,p),1,max))})
    } else{
      gaga <- apply(ans, 2, function(x){median(apply(matrix(x,n.err,p),1,max))})
    }
  }

  return(list(    gaga = gaga,
              messages = messages))
}
