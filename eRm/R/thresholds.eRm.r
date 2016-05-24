thresholds.eRm <- function(object)                # uses matrix approach
{
#Computation of threshold parameters for polytomous models
#object of class "eRm" (but not "dRm")

  if(object$model %in% c("LLTM", "RM")) stop("Threshold parameters are computed only for polytomous models!")
  if(object$model %in% c("LRSM", "LPCM")) {
    mpoints <- object$mpoints
    ngroups <- object$ngroups
    vecrep  <- mpoints * ngroups                      
  } else {
    mpoints <- 1L
    ngroups <- 1L
    vecrep  <- 1L
  }
  
  betapar    <- object$betapar
  indmt      <- apply(object$X, 2L, max, na.rm = TRUE)      # number of categories per item
  mt_vek1    <- sequence(indmt[1L:(length(indmt)/mpoints)]) # 1 block of beta-items
  mt_vek     <- rep(mt_vek1, vecrep)
  sq         <- ifelse(mt_vek > 1, -1, 0)
  d1         <- diag(sq[-1L])
  k          <- length(betapar)
  d2         <- diag(k)
  d2[-k,-1L] <- d2[-k, -1L] + d1
  T_mat      <- t(d2)                           # MM 2010-02-20
  threshpar  <- -as.vector(T_mat %*% betapar)   # vector with threshold parameters - fix: MM 2010-02-20
  
  names(threshpar) <- paste("thresh", names(betapar))
  
  vc.beta   <- object$W %*% solve(object$hessian) %*% t(object$W) # VC matrix beta's
  se.thresh <- sqrt(diag( T_mat %*% vc.beta %*% t(T_mat) ))       # standard errors of thresholds - fix: MM 2010-02-20
  names(se.thresh) <- names(threshpar)

  blocks  <- rep(1L:vecrep, each = length(mt_vek1))
  thblock <- split(threshpar, blocks)                 #block of threshholds (as in design matrix)
  indmt1  <- indmt[1L:(length(indmt)/mpoints)]
  indvec  <- rep(1L:length(indmt1), indmt1)

  threshtab.l <- lapply(thblock, function(x) {                     #list of table-blocks
                     location  <- tapply(x,indvec,mean)             #location parameters
                     thresh.l  <- split(x, indvec)
                     threshmat <- t(as.data.frame(lapply(thresh.l, function(i_th){
                                    c(i_th, rep(NA, length.out=max(mt_vek)-length(i_th)))
                                  })))
                     colnames(threshmat) <- paste("Threshold", 1:dim(threshmat)[2])
                     parmat <- cbind("Location" = location, threshmat)
  }) 
  
  #determine item names for block-table
  cnames   <- colnames(object$X)
  ind.it   <- rep(1L:mpoints, each = length(cnames)/mpoints)           #item label index
  itnames1 <- as.vector(unlist(tapply(cnames, ind.it, function(x){ rep(x, ngroups) }))) 
  rep.ind  <- unlist(lapply(threshtab.l, nrow))
  sp.ind   <- rep(1L:length(rep.ind), rep.ind)

  names.l <- split(itnames1, sp.ind)                   #names as list
  for(i in seq_along(threshtab.l)) rownames(threshtab.l[[i]]) <- names.l[[i]]              #name the items

  result <- list("threshpar"   = threshpar,
                 "se.thresh"   = se.thresh,
                 "threshtable" = threshtab.l)
  class(result) <- "threshold"

  return(result)

}
