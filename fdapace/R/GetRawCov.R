GetRawCov <- function(y,t,obsGridnew, mu, dataType, error){
#  obtain raw covariance
#  Input y :       1*n cell array of the observed repeated measurements from n subjects
#  Input t :       1*n cell array of the observed time points from n subjects
#  Input obsGridnew:  1*m vector of time points correspond to mu
#  Input mu:       1*m vector of fitted mean functions from Step I, corresponding to
#                 pooled unique time points from t
#  Input dataType: Output of IsRegular()
#  Input error:    TRUE with measurement error assumption
#                  FALSE without measurement error assumption
# 
#  Output res: a list that contains tPairs, cxxn, indx,win and cyy
#     tPairs:  N * 2  vector denotes the  pairs of time points for subject 
#                 concatenating as two vectors 
#                if error = 1, all (t_ij, t_ij) will be removed  
#       cxxn:    1 * N vector of raw covariance corresponding to tPairs      
#       indx:    1 * N vector of indices for each subject
#        win:    1 * N weight matrix for the 2-D smoother for covariance function
#        cyy:    1 * M vector of raw covariance corresponding to all pairs of time points,
#                i.e., it is the same as cxxn if error = 0
#       diag:    if error == TRUE: 2-column matrix recording raw covariance along the diagonal direction (col 2)
#                and the corresponding observed time points (col 1)   
#                if error == FALSE: NULL

  ncohort <- length(y);
  obsGrid <- sort(unique(unlist(t)))
  mu <- MapX1D(x = obsGridnew, y = mu, newx = obsGrid);
  count <- NULL
  indx = NULL 
  diag = NULL

  if(dataType == 'Sparse'){
  
    Ys = lapply(X = y, FUN=pracma::meshgrid) #pracma
    Xs = lapply(X = t, FUN=pracma::meshgrid) #pracma

    # vectorise the grids for y & t
    xx1 = unlist(do.call(rbind, lapply(Xs, '[', 'X')) )
    xx2 = unlist(do.call(rbind, lapply(Xs, '[', 'Y')) ) 
    yy2 = unlist(do.call(rbind, lapply(Ys, '[', 'Y')) )
    yy1 = unlist(do.call(rbind, lapply(Ys, '[', 'X')) )
    
    # get id1/2 such that xx1/2 = q(id1/2), where q = unique(xx1/2)
    # id1 = apply(X= sapply(X=xx1, FUN='==',  ...=sort(unique(xx1)) ),MARGIN=2, FUN=which)
    # id2 = apply(X= sapply(X=xx2, FUN='==',  ...=sort(unique(xx2)) ),MARGIN=2, FUN=which)
    # This is more stable and faster.
    id1 = uniqueM(xx1)
    id2 = uniqueM(xx2)
    cyy = ( yy1 - mu[ id1]) * (yy2 - mu[id2] )
    
    # index for subject i
    # indx = unlist(sapply( 1:length(y), function(x) rep(x,  (unlist(lapply(length, X= y))[x])^2) ))
    # This is more stable and faster.
      indx = rep( 1:length(y), times =  unlist(lapply(y, length))^2)

    tPairs = matrix( c(xx1, xx2), nrow=length(xx1), ncol=2);

    if(error){
      tneq = which(xx1 != xx2)
      teq = which(xx1 == xx2)
      indx = indx[tneq];
      diag = matrix(c(tPairs[teq,1], cyy[teq]), ncol = 2)
      tPairs = tPairs[tneq,];
      cxxn = cyy[tneq];     
    }else{
      cxxn = cyy;     
    }

    # win = pracma::ones(1, length(cxxn));
    # count = GetCount(tPairs)...

  }else if(dataType == 'Dense'){
    
    yy = t(matrix(unlist(y), length(y[[1]]), ncohort))
    MU = t(matrix( rep(mu, times=length(y)), ncol=length(y)))
    t1 = t[[1]]
    yy = yy - MU;
    cyy = t(yy) %*% yy / ncohort
    cyy = as.vector(t(cyy))
    cxxn = cyy;
    xxyy = pracma::meshgrid(t1); # pracma

    tPairs =  (matrix( c(c(xxyy$X), c(xxyy$Y)), ncol = 2))

    if(error){
      tneq = which(tPairs[,1] != tPairs[,2])
      teq = which(tPairs[,1] == tPairs[,2])
      diag = matrix(c(tPairs[teq,1], cyy[teq]), ncol = 2)
      tPairs = tPairs[tneq,];
      cxxn = cyy[tneq];     
    }else{
      cxxn = cyy;     
    }

   # win = pracma::ones(1, length(cxxn));
  }else if(dataType == 'RegularWithMV'){
    stop("This is not implemented yet. Contact Pantelis!")
  }else {
    stop("Invalid 'dataType' argument type")
  } 
    
  result <- list( 'tPairs' = tPairs, 'cxxn' = cxxn, 'indx' = indx, # 'win' = win,
   'cyy' = cyy, 'diag' = diag, 'count' = count, 'error' = error, 'dataType' = dataType);  
 
  class(result) <- "RawCov"
  return(result)
}

 uniqueM <- function(x){ 
   g = sort(unique(x))
   id1 = rep(0, length(x));
   for( i in 1:length(g)){
    id1[which(x == g[i])] = i;
   }
   return(id1)
 }


