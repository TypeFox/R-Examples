GetRawCrCovFuncFunc <- function(Ly1, Lt1 = NULL, Ymu1, Ly2, Lt2 = NULL, Ymu2){

  # If Ly1 and Ly2 are both matrices and Lt1 and Lt2 are both null then assume DENSE
  if ( is.matrix(Ly1) && is.matrix(Ly2) && is.null(Lt1) && is.null(Lt2) ){
     if( dim(Ly1)[1] != dim(Ly2)[1] ){
       stop("Ly1 and Ly2 are not compatible")
     } 
     KK = cov(Ly1,Ly2,use="pairwise.complete.obs")
     return( list(rawCCov =  1 * (KK), tpairn = NULL) )    
  }
  # otherwise assume SPARSE
  if ( ! CheckEqualLengths(Lt1,Ly1)){
    stop("Lt1 and Ly1 are not compatible")
  }
  if ( ! CheckEqualLengths(Lt2,Ly2) ){
    stop("Lt2 and Ly2 are not compatible")
  }    
  ulLt1 = unlist(Lt1)
  ulLt2 = unlist(Lt2)
  if ( ! CheckEqualLengths(unique(ulLt1), Ymu1)){
    stop("Lt1 and Ymu1 are not compatible")
  }       
  if ( ! CheckEqualLengths(unique(ulLt2), Ymu2)){
    stop("Lt2 and Ymu2 are not compatible")
  }       

  # Centre both lists according to their means    
  muY1 <- approxfun(x= sort(unique(ulLt1)), y = Ymu1)
  Ly1c <- lapply(1:length(Ly1), function(i) Ly1[[i]]- muY1( Lt1[[i]]) )
  muY2 <- approxfun(x= sort(unique(ulLt2)), y = Ymu2)
  Ly2c <- lapply(1:length(Ly2), function(i) Ly2[[i]]- muY2( Lt2[[i]]) )

 # I keep this so the functional code below is understandable in iterative form
 # tPairs1 <- c()
 # tPairs2 <- c()
 # cyy <- c()
 # for (i in 1:length(Ly1c)){
 #   q = length(Ly1c[[i]])
 #   p = length(Ly2c[[i]])
 #   cyy <- c( cyy, rep(x= Ly1c[[i]],each=p) * rep(x= Ly2c[[i]],times=q)  )  
 #   tPairs1 <- c( tPairs1,  rep(Lt1[[i]],each=p))
 #   tPairs2 <- c( tPairs2,  rep(Lt2[[i]],times=q))
 # }

  cyy <- as.vector(unlist (mapply(FUN = function(v1, v2){ return( rep(v1,times = length(v2)) * rep(v2, each = length(v1) )) }, v1= Ly2c, v2 = Ly1c)) )
  tPairs2 <- as.vector(unlist (mapply(FUN = function(v1, v2){ return(rep(v1,times = length(v2))) }, v1= Lt2, v2 = Ly1c)) )
  tPairs1 <- as.vector(unlist (mapply(FUN = function(v1, v2){ return(rep(v1, each = length(v2))) }, v1= Lt1, v2 = Ly2c)) )

  RCC <-list()
  RCC$rawCCov = cyy
  RCC$tpairn = cbind(tPairs1, tPairs2)
  RCC$IDs = rep( 1:length(Lt1), times = sapply(Lt1, length) * sapply(Lt2, length)  )
  return(RCC)
}

CheckEqualLengths <- function(x1,x2){ return( all.equal( sapply(x1, length) , sapply(x2, length) )) } 
