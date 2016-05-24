#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# SSVD.iter.thresh : FIT-SSVD method to obtain SVD                             #
#    Original code provided by Dr. Dan Yang, modified only in commentary and   #
#    treatment of warnings as messages.                                        #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x            : matrix of observed data                                     #
#   u.old        : original matrix of left singular vectors                    #
#   v.old        : original matrix of right singular vectors                   #
#   gamma.u      : constant used in simple threshold method                    #
#   gamma.v      : constant used in simple threshold method                    #
#   dothres      : type of threshold (soft/hard)                               #
#   num.sig      : number of elements in threshold                             #
#   tol          : tolerance level for distance between subspaces              #
#   n.step       : maximum number of iterations allowed                        #
#   sigma        : standard deviation                                          #
#   methodology  : use Algorithm 3 to obtain thresholds                        #
#   n.err        : number of bootstraps in Algorithm 3                         #
#   error.full   : error estimate method - only full available                 #
#   error.median : use median in Algorithm 3                                   #
# Outputs                                                                      #
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

SSVD.iter.thresh <- function(x, 
                             u.old,  
                             v.old,  
                             gamma.u=sqrt(2),  
                             gamma.v=sqrt(2), 
                             dothres="hard", 
                             num.sig=dim(u.old)[2],  
                             tol=1e-8,  
                             n.step=100,  
                             sigma=NA,  
                             methodology=TRUE, 
                             n.err=100,  
                             error.full=TRUE,  
                             error.median=TRUE){

  messages <- list()

  x <- as.matrix(x)
  p <- dim(x)[1]
  q <- dim(x)[2]
  pq <- max(p,q)

  # estimate standard deviation
  if(is.na(sigma)){
    sigma.hat  <-  mad(as.vector(x))
  } else {
    sigma.hat <- sigma
  }

  # rescale X by standard deviation
  x.scaled <- x/sigma.hat

  # identify thresholding function
  if(dothres == "hard"){
    thresh <- hard.thresh
  } else if(dothres == "soft"){
    thresh <- soft.thresh
  } else {
    messages <- c(messages, 
      "SSVD.iter.thresh: argument dothres not recognized! Hard-thresholding used as default")
    thresh <- hard.thresh
  }

  # identify error function
  if(error.full == TRUE){
    error.est <- error.est.full
  } else {
#    error.est <- error.est.row.col
    error.est <- error.est.full  # error.est.row.col was not provided in original code set
  }

  # initializations
  dist.u <- 1
  dist.v <- 1
  i.step <- 0
  u.cur <- u.old
  v.cur <- v.old
  threshold.v <- rep(NA, num.sig)

  while(i.step <= n.step & max(dist.u,dist.v) > tol ){

    u.old <- u.cur

    #----------------------------------------------------------------------#
    # Right-to-left Multiplication                                         #
    #----------------------------------------------------------------------#
    u.cur <- x.scaled %*% v.cur

    #----------------------------------------------------------------------#
    # Left Thresholding                                                    #
    #----------------------------------------------------------------------#
    if(methodology == FALSE){
      u.cur <- thresh(u.cur, gamma.u*sqrt(log(pq)))
      threshold.u <- gamma.u*sqrt(log(pq)) * rep(1, num.sig)
    } else {
      threshold.u <- error.est(x.scaled, u.old, v.cur, n.err=n.err, 
                               error.median=error.median)
      messages <- c(messages, threshold.u$messages)
      threshold.u <- threshold.u$gaga
      u.cur <- thresh(u.cur, threshold.u)
    }

    #----------------------------------------------------------------------#
    # Left Orthonormalization with QR Decomposition                        #
    #----------------------------------------------------------------------#
    if(all(u.cur==0)){
      messages <- c(messages, "SSVD.iter.thresh: Overthresh: all zero!")
      u.cur <- u.old
      dist.u <- 0
      dist.v <- 0
      break
    } else {
      u.cur <- qr.Q(qr(u.cur))
      dist.u <- subsp.dist.orth(u.cur,u.old)
    }

    v.old <- v.cur

    #----------------------------------------------------------------------#
    # Left-to-Right Multiplication                                         #
    #----------------------------------------------------------------------#
    v.cur <- t(x.scaled) %*% u.cur

    #----------------------------------------------------------------------#
    # Right Thresholding                                                   #
    #----------------------------------------------------------------------#
    if(methodology == FALSE){
      v.cur <- thresh(v.cur, gamma.v*sqrt(log(pq)))
      threshold.v <- gamma.v*sqrt(log(pq)) * rep(1, num.sig)
    } else {
      threshold.v <- error.est(t(x.scaled), v.old, u.cur, 
        n.err=n.err, error.median=error.median)
      messages <- c(messages, threshold.v$messages)
      threshold.v <- threshold.v$gaga
      v.cur <- thresh(v.cur, threshold.v)
    }

    #----------------------------------------------------------------------#
    # Right Orthonormalization with QR Decomposition                       #
    #----------------------------------------------------------------------#
    if(all(v.cur==0)){
      messages <- c(messages, "SSVD.iter.thresh: Overthresh: all zero!")
      v.cur <- v.old
      dist.u <- 0
      dist.v <- 0
      break
    } else {
      v.cur <- qr.Q(qr(v.cur))
      dist.v <- subsp.dist.orth(v.cur,v.old)
    }

    i.step <- i.step + 1
  }

  k.u <- apply(u.cur!=0, 2, sum)
  k.v <- apply(v.cur!=0, 2, sum)
  d.cur <- diag(t(u.cur) %*% x %*% v.cur)
  ans <- consistent.signs(u.cur,d.cur,v.cur)

  if(i.step==n.step){
    messages <- c(messages, "maximum iterations used; increase nstep")
  }

  list(        "u" = ans$u, 
               "v" = ans$v, 
               "d" = ans$d, 
             "k.u" = k.u, 
             "k.v" = k.v, 
           "nstep" = i.step - 1,
       "sigma.hat" = sigma.hat,
           "thr.u" = threshold.u,
           "thr.v" = threshold.v,
            "dist" = max(dist.u,dist.v),
        "messages" = messages)
}
