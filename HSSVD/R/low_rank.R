#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# low_rank : main function of HSSVD package                                    #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x      : Matrix of observed data                                           #
#   ranks  : Optional : integer vector (dim=3) specifying rank estimates for   #
#            each SVD step of the framework (Steps 2, 3 & 4)                   #
# ---------------------------------------------------------------------------- #
# if rank estimates not provided by user                                       #
# ---------------------------------------------------------------------------- #
#   sparse : logical vector (dim=3) indicating matrix at each step is sparse   #
#            Default: (FALSE,FALSE,FALSE)                                      #
#   est    : character vector (dim=3) specifying rank estimation method,       #
#            only options are 'wold' or 'gb'.  Default = 'wold'                #
#   add    : integer vector (dim=3) 0 <= add[i] <= 2 to be added to rank       #
#            estimate. Default: (1,0,0)                                        #
#   ...    : additional inputs to be passed to FIT.SSVD                        #
# Outputs                                                                      #
#   call        : record of call to low_rank                                   #
#   rescale     : SVD of Approximation matrix U (Step 2)                       #
#   result_mean : SVD of Approximation matrix Y_tilde (Step 3)                 #
#   result_var  : SVD of Approximation matrix Z_tilde (Step 4)                 #
#   ranks       : vector of estimated ranks for (U, Y_tilde, Z_tilde)          #
#   bgmean      : mean of null cluster (Step 5)                                #
#   bgstd       : standard deviation of null cluster (Step 5)                  #
#   back        : logicial matrix; background cluster (TRUE) (Step 5)          #
#   mean_app    : mean approximation (Step 6)                                  #
#   std_app     : standard deviation approximation (Step 6)                    #
#   messages    : Informative messages provided during calculation             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
low_rank <- function(x,
                     ranks=c(rep(NA, 3)),
                     sparse=rep(FALSE,3),
                     est=rep("wold",3),
                     add=c(1,0,0), ...){


  call <- match.call()

  #--------------------------------------------------------------------------#
  # Verify and initialize input                                              #
  #--------------------------------------------------------------------------#
  est <- tolower(est)

  messages <- list()

  #--------------------------------------------------------------------------#
  # All ranks are given or estimated - not a mix and match                   #
  #--------------------------------------------------------------------------#
  if( (sum(is.na(ranks)) != 0) && (sum(is.na(ranks)) != 3 ) ){
    stop("If provided, ranks must be provided for all steps.")
  }

  #--------------------------------------------------------------------------#
  # If given, all ranks must be greater than 0                               #
  #--------------------------------------------------------------------------#
  if(!any(is.na(ranks))){
    if(any(ranks<=0)){
      stop("ranks must be > 0")
    }
  }

  #--------------------------------------------------------------------------#
  # ranks or (est, sprase and add) must be provided for each step            #
  #--------------------------------------------------------------------------#
  if(length(ranks) != 3 || length(sparse) != 3 || 
     length(est)   != 3 || length(add)    != 3){
    stop("Incorrect dimensions for input variables. Must be dim=3")
  }

  #--------------------------------------------------------------------------#
  # verify the estimation type provided - must be 'gb' or 'wold'             #
  #--------------------------------------------------------------------------#
  if(any(!(est %in% c("gb","wold")))){
    stop("Verify est input variable.")
  }

  #--------------------------------------------------------------------------#
  # all add must be between 0 and 2                                          #
  #--------------------------------------------------------------------------#
  if(any(add > 2) || any(add < 0)){
    stop("add must be 0 <= add <= 2")
  }

  arg.list <- list(...)
  if(length(arg.list) > 0){
    it <- names(arg.list) %in% c('dothres', 'n.step', 'n.err')
    arg.list[!it] <- NULL
  }

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 1::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Standardize x                                                            #
  #--------------------------------------------------------------------------#
  X.std <- (x - median(x))/mad(x)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 2::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Quadratic Rescaling                                                      #
  #--------------------------------------------------------------------------#
  X.scaled <- X.std^2 - 1

  #--------------------------------------------------------------------------#
  # Apply SSVD on X.scaled                                                   #
  #--------------------------------------------------------------------------#

  tmp <- onestep(X.scaled, sparse[1], est[1], ranks[1], add[1], arg.list)
  X.scaled.SSVD <- tmp$result
  X.scaled.rank <- tmp$rank
  messages <- c(messages, tmp$messages)

  #--------------------------------------------------------------------------#
  # Proportion of variance explained by the top i rank sparse matrix, R_i    #
  # If the new rank estimated as R_r > 0.95 1 <= r <= (k+1) is less than     #
  # original rank, recalculate.                                              #
  #--------------------------------------------------------------------------#
  rank_new <- sum(cumsum(X.scaled.SSVD$d^2) < 0.95*sum(X.scaled.SSVD$d^2))+1
  if(rank_new < X.scaled.rank){
    X.scaled.rank <- rank_new
    eigenX <- svd(X.scaled)
    X.scaled.SSVD <- layer.find(x=X.scaled,
                     u0=matrix(eigenX$u[,1:X.scaled.rank],ncol=X.scaled.rank),
                     v0=matrix(eigenX$v[,1:X.scaled.rank],ncol=X.scaled.rank),
                     arg.list)
    messages <- c(messages, X.scaled.SSVD$messages)
    X.scaled.SSVD$messages <- NULL
  }

  est.U <- X.scaled.SSVD$u %*% 
           diag(X.scaled.SSVD$d,ncol=X.scaled.rank,nrow=X.scaled.rank) %*% 
           t(X.scaled.SSVD$v)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 3::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Define Y                                                                 #
  #--------------------------------------------------------------------------#
  est.var <- sqrt(est.U + 1 - min(1.001*(est.U+1),0))
  Y <- X.std/est.var

  #--------------------------------------------------------------------------#
  # Apply SSVD on Y                                                          #
  #--------------------------------------------------------------------------#
  tmp <- onestep(Y, sparse[2], est[2], ranks[2], add[2], arg.list)
  Y.SSVD <- tmp$result
  Y.rank <- tmp$rank
  messages <- c(messages, tmp$messages)

  #--------------------------------------------------------------------------#
  # Proportion of variance explained by the top i rank sparse matrix, R_i    #
  # If the new rank estimated as R_r > 0.95 1 <= r <= (k+1) is less than     #
  # original rank, recalculate.                                              #
  #--------------------------------------------------------------------------#
  rank_new <- sum(cumsum(Y.SSVD$d^2) < 0.95*sum(Y.SSVD$d^2)) + 1
  if(rank_new < Y.rank){
   Y.rank <- rank_new
   eigenY <- svd(Y)
   Y.SSVD <- layer.find(x=Y,
               u0=matrix(eigenY$u[,1:Y.rank],ncol=Y.rank),
               v0=matrix(eigenY$v[,1:Y.rank],ncol=Y.rank),
               arg.list)
    messages <- c(messages, Y.SSVD$messages)
    Y.SSVD$messages <- NULL
  }

  est.Y <- Y.SSVD$u %*% diag(Y.SSVD$d,ncol=Y.rank,nrow=Y.rank) %*% t(Y.SSVD$v)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 4::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Define Z                                                                 #
  #--------------------------------------------------------------------------#
  Z_origin <- (X.std - est.Y*est.var)
  Z_origin[Z_origin==0] <- min(abs(Z_origin[abs(Z_origin > 0)]))
  Z <- log(Z_origin^2) - mean(log(Z_origin^2))

  #--------------------------------------------------------------------------#
  # Apply SSVD on Z                                                          #
  #--------------------------------------------------------------------------#
  tmp <- onestep(Z, sparse[3], est[3], ranks[3], add[3], arg.list)
  Z.SSVD <- tmp$result
  Z.rank <- tmp$rank
  messages <- c(messages, tmp$messages)

  #--------------------------------------------------------------------------#
  # Proportion of variance explained by the top i rank sparse matrix, R_i    #
  # If the new rank estimated as R_r > 0.95 1 <= r <= (k+1) is less than     #
  # original rank, recalculate.                                              #
  #--------------------------------------------------------------------------#
  rank_new <- sum(cumsum(Z.SSVD$d^2) < 0.95*sum(Z.SSVD$d^2)) + 1
  if(rank_new < Z.rank){
    Z.rank <- rank_new
    eigenZ <- svd(Z)
    Z.SSVD <- layer.find(x=Z,
                u0=matrix(eigenZ$u[,1:Z.rank],ncol=Z.rank),
                v0=matrix(eigenZ$v[,1:Z.rank],ncol=Z.rank),
                arg.list)
    messages <- c(messages, Z.SSVD$messages)
    Z.SSVD$messages <- NULL
  }

  est.Z <- Z.SSVD$u %*% diag(Z.SSVD$d,ncol=Z.rank,nrow=Z.rank) %*% t(Z.SSVD$v)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 5::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Define matrix of indicators                                              #
  #--------------------------------------------------------------------------#
  P <- (est.Y == 0) & (est.Z == 0)
  #--------------------------------------------------------------------------#
  # estimate b_hat and rho_hat                                               #
  #--------------------------------------------------------------------------#
  bgmean <- mean(x[P])
  bgstd <- sd(x[P])

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::STEP 6::::::::::::::::::::::::::::::::::::#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #--------------------------------------------------------------------------#
  # Define matrix of indicators for est.Y==0                                 #
  #--------------------------------------------------------------------------#
  P1 <- est.Y == 0
  #--------------------------------------------------------------------------#
  # calculate mean approximation                                             #
  #--------------------------------------------------------------------------#
  my.mean <- mad(x)*est.Y*est.var + median(x)*(1-P1) + bgmean*P1

  #--------------------------------------------------------------------------#
  # Define matrix of indicators for est.Z==0                                 #
  #--------------------------------------------------------------------------#
  P2 <- est.Z == 0
  #--------------------------------------------------------------------------#
  # calculate standard deviation approximation                               #
  #--------------------------------------------------------------------------#
  my.var <- (bgstd^2*P2 + mad(x)^2*(1-P2))*exp(est.Z)
  my.sd <- sqrt(my.var)

  return(list(        call = call,
                   rescale = X.scaled.SSVD,
               result_mean = Y.SSVD,
                result_var = Z.SSVD,
                     ranks = c(X.scaled.rank, Y.rank, Z.rank),
                    bgmean = bgmean,
                     bgstd = bgstd,
                      back = P,
                  mean_app = my.mean,
                   std_app = my.sd,
                  messages = messages))

}
