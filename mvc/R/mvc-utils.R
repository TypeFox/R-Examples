#' Check views for consistency
#' Views must have exactly the same row names
#' @param view1 View 1
#' @param view2 View 2
#' @return Message and stop as appropriate

checkViews <- function(view1, view2) {
  if (!is.data.frame(view1) | !is.data.frame(view2)) {
    stop("Preparing data failed: Views must be data frames")
  }
  if ( sum(sort(row.names(view1))==sort(row.names(view2)) ) != length(row.names(view1)) ) {
    stop("Preparing data failed: Row names differ")
  }
  if ( ! identical( row.names(view1), row.names(view2)) ) {
    stop("Preparing data failed: Order of row names differs between views")
  }

}


#' Counts unique values in both views
#' Stops on any non-numeric values
#' @param view1 View 1
#' @param view2 View 2
#' @return list containing unique values for each view

viewsClasses <- function(view1, view2) {
  nrColsView1 <- dim(view1)[2]
  nrColsView2 <- dim(view2)[2]
  if (sum(apply(view1,2,is.numeric)) < nrColsView1) {
    stop("view1 contains non-numeric columns")
  }
  if (sum(apply(view2,2,is.numeric)) < nrColsView2) {
    stop("view2 contains non-numeric columns")
  }
  uniqueValsView1 <- sort(unique( unlist(as.list (apply(view1,2,unique),use.names=F) ) ))
  uniqueValsView2 <- sort(unique( unlist(as.list (apply(view2,2,unique),use.names=F) ) ))
  list(view1=uniqueValsView1, view2=uniqueValsView2)
}



# # # Multi-Spherical k-Means # # #


#' Euclidean length of vector
#' @param x vector
#' @return length of x

vectorLength <- function(x) {
  as.vector(sqrt(x %*% x))
}



#' Unit length for vector
#' @param x vector
#' @return x converted to unit length

UL <- function(x) {
  x = x/vectorLength(x) # trick to scale to unit length
  x[is.nan(x)] = 0
  x
}


#' Unit length of all vectors row-wise
#' @param X matrix
#' @return X row-wise converted to unit length

rowWUL <- function(X) {
  # X = t ( apply ( X, 1, function(x) x/vectorLength(x) ) ) # trick to scale to unit length
  for (i in 1:NROW(X)) {
    #tryCatch ( X[i,]=X[i,]/vectorLength(X[i,]), error=function(e) print(X[i,]) ) 
    X[i,]=X[i,]/vectorLength(X[i,])
  }
  X[is.nan(X)] = 0
  X
}


#' Objective Function (sum of cosines)
#' @param X data matrix (row-wise vectors in unit length).
#' @param C concept vectors as matrix (row-wise in unit length).
#' @param CIdx vector of length NROW(X) with natural numbers 1..k, indicating cluster for each data vector.
#' @return sum of cosine-similarities.
#' @examples { 
#'   X=structure(c(0.707, 0.707, 0.707, 0.707), .Dim = c(2L, 2L))
#'   C=structure(c(1, 0, 0, 1), .Dim = c(2L, 2L))
#'   CIdx=c(2, 1)
#'   oFSkm(X,C,CIdx) # 1.414
#' }

oFSkm <- function(X,C,CIdx) {
  CIdxs = unique(sort(CIdx))
  k = length(CIdxs)
  of = vector( length = k, mode = 'numeric' ) # mode 'numeric' inits to 0
  for ( j in 1:k ) {
    subSpaceJ = X [ CIdx == CIdxs[j], , drop=F ]
    #if (dim(subSpaceJ)[1] > 0) {
      for ( i in 1:NROW(subSpaceJ) ) {
        of[j] = of[j] + abs( unlist(subSpaceJ[i,]) %*% unlist(C[j,]) )
      }
    #}
  }
  sum(of)
}


#' Calculate concept vectors for Spherical k-Means as unit length sum of vectors of the k clusters.
#' @param X data matrix (row-wise in unit length).
#' @param CIdx vector of length NROW(X) with natural numbers 1..k, indicating cluster for each data vector.
#' @param doOutput whether progress bar indicators should be output
#' @return concept vectors as matrix (row-wise in unit length).
#' @examples {
#'   X=structure(c(1, 1, -1, 0, 1, 0, -1, -1), .Dim = c(4L, 2L))
#'   CIdx=c(1, 1, 2, 2)
#'   C=conceptVectorsSkm(X,CIdx)
#'   dput(C) 
#'   # structure(c(0.894427190999916, -0.447213595499958, 
#'   # 0.447213595499958, -0.894427190999916), .Dim = c(2L, 2L))
#' }

conceptVectorsSkm <- function(X,CIdx,doOutput=F) {
  CIdxs = unique(sort(CIdx))
  k = length(CIdxs)
  C = matrix( rep( 0, k*NCOL(X) ), k, NCOL(X) )
  if (doOutput) pb <- txtProgressBar(min = 0, max = k, style = 3)
  for (i in 1:k) {
    ithConceptSum = colSums( X [ CIdx==CIdxs[i], , drop=F ] )
    C[i,] = UL(unlist(ithConceptSum))
    if (doOutput) setTxtProgressBar(pb, i)
  }
  if (doOutput) close(pb)
  C
}


#' Calculate partitions (concept indices) by assigning each vector to the closest concept vector.
#' @param X data matrix (row-wise in unit length).
#' @param C matrix with k rows, indicating concept vectors (row-wise in unit length).
#' @param doOutput whether progress bar indicators should be output
#' @return concept indices as vector.
#' @examples { 
#'   X=structure(c(1, 1, -1, 0, 1, 0, -1, -1), .Dim = c(4L, 2L))
#'   C=structure(c(0.894427190999916, -0.447213595499958, 
#'   0.447213595499958, -0.894427190999916), .Dim = c(2L, 2L))
#'   CIdx=conceptIndicesSkm(X,C)
#'   dput(CIdx) 
#'   # c(1, 1, 2, 2)
#' }

conceptIndicesSkm <- function(X,C,doOutput=F) {
  k = NROW(C)
  n = NROW(X)
  CIdx = rep( 0, n )
  if (doOutput) pb <- txtProgressBar(min = 0, max = n, style = 3)
  for ( i in 1:n ) {
    maxDist = -Inf # Init to -infty
    for ( j in 1:k ) {
      actualDist = abs(as.vector(unlist(X[i,]) %*% unlist(C[j,])))
      if (actualDist > maxDist) {
        maxDist = actualDist
        CIdx[i] = j
      }
    }
    if (doOutput) setTxtProgressBar(pb, i)
  }
  if (doOutput) close(pb)
  CIdx
}


#' Calculate means per Cluster and view for Spherical k-Means by using a consensus approach.
#' @param view1 data matrix (row-wise in unit length).
#' @param view2 data matrix (row-wise in unit length).
#' @param view1Idx vector of length NROW(view1) with natural numbers 1..k, indicating cluster for each data vector of view1.
#' @param view2Idx vector of length NROW(view1) with natural numbers 1..k, indicating cluster for each data vector of view2.
#' @return cluster means as matrices per view (row-wise in unit length).
#' @examples {
#'   view1 = structure(c(1, 1, -1, 0, 1, 0, -1, -1), .Dim = c(4L, 2L))
#'   view2 = structure(c(1, 1, -1, 0, 1, 0, -1, 0), .Dim = c(4L, 2L))
#'   view1Idx = c(2, 2, 1, 1)
#'   view2Idx = c(2, 1, 1, 1)
#'   mPerClV=consensusMeansPerClVSkm(view1,view2,view1Idx,view2Idx)
#'   dput(mPerClV) 
#' }

consensusMeansPerClVSkm <- function(view1,view2,view1Idx,view2Idx) {

# structure(list(view1 = structure(c(-0.447213595499958, 0.707106781186547, -0.894427190999916, 0.707106781186547), .Dim = c(2L, 2L)), view2 = structure(c(-0.707106781186547, 0.707106781186547, -0.707106781186547, 0.707106781186547), .Dim = c(2L, 2L))), .Names = c("view1", "view2"))

  #n = length(view1Idx)

  intersectOfIdxs = intersect(sort(unique(view1Idx)),sort(unique(view2Idx)))
  k = length(intersectOfIdxs)

  #if (sum(sort(unique(view1Idx)) != sort(unique(view2Idx)))>0) {
  #  print(unique(view1Idx))
  #  print(unique(view2Idx))
  #  stop("Different clusters.")
  #}

  view1mMatrix = matrix(rep(0,k*NCOL(view1)),k,NCOL(view1),byrow=T)
  view2mMatrix = matrix(rep(0,k*NCOL(view2)),k,NCOL(view2),byrow=T)

  for ( j in 1:k ) {
    sharedInJ = ((view1Idx==intersectOfIdxs[j]) & (view2Idx==intersectOfIdxs[j]))
    mPerClV1 = UL( unlist( colSums( view1 [sharedInJ, ,drop=F] ) ) )
    view1mMatrix[j,]=mPerClV1
    mPerClV2 = UL( unlist( colSums( view2 [sharedInJ, ,drop=F] ) ) )
    view2mMatrix[j,]=mPerClV2
  }

  list("view1"=view1mMatrix, "view2"=view2mMatrix)
}


#' Assign final indices to means that have the smallest angle.
#' @param view1 data matrices (row-wise in unit length).
#' @param view2 data matrices (row-wise in unit length).
#' @param mPerClV list of means per Cluster and View.
#' @return vector of indices for each data point.
#' @examples \dontrun{ 
#'   view1 = structure(c(1, 1, -1, 0, 1, 0, -1, -1), .Dim = c(4L, 2L))
#'   view2 = structure(c(1, 1, -1, 0, 1, 0, -1, 0), .Dim = c(4L, 2L))
#'   finIdx = assignFinIdxPerClSkm(view1,view2,mPerClV)
#'   dput(finIdx) 
#'   # c(2, 2, 1, 1)
#' }

assignFinIdxPerClSkm <- function(view1,view2,mPerClV) {

  view1 = matrix ( unlist(view1), dim(view1)[1], dim(view1)[2] )
  view2 = matrix ( unlist(view2), dim(view2)[1], dim(view2)[2] )
  view1 = rowWUL(view1)
  view2 = rowWUL(view2)
  mPerClV1 = rowWUL(mPerClV[["view1"]])
  mPerClV2 = rowWUL(mPerClV[["view2"]])

  k = NROW(mPerClV1)
  if (NROW(mPerClV2) != k) stop("Different clusters.")

  n = NROW(view1)
  finalCIdx = rep(0,n)
  sphericMin = rep(Inf,n) # init to +infty
  for ( j in 1:k ) {
    sphericVal = as.vector( acos( round( view1 %*%  mPerClV1[j,] , digits=5 ) ) + acos( round(  view2 %*% mPerClV2[j,] , digits=5 ) ) )
    minPat = sphericVal < sphericMin
    sphericMin[minPat] = sphericVal[minPat]
    finalCIdx[minPat] = j
  }
  finalCIdx
}



# # # Mixture of Categoricals EM # # #


#' Calculate Bernoulli likelihood
#' @param x a binary event (vector)
#' @param prob the Bernoulli probability (vector)
#' @return Bernoulli likelihood

dbern <- function(x,prob) {
  prob^x * (1-prob)^(1-x)
}

#' Calculate categorical likelihood
#' @param x a categorical event vector
#' @param prob the categorical probability matrix (rows along events, cols along event values)
#' @return categorical likelihood
#' @examples {
#'   dcat(c(1,2,1),matrix(c(.9,.8,.9,.1,.2,.1),3,2))
#' }

dcat <- function(x,prob) {
    sapply(seq_along(x), function(idx) prob[idx,x[idx]])
}


#' Calculate Bernoulli likelihood row-wise for binary events
#' @param X a matrix of binary events (row-wise)
#' @param prob the Bernoulli probability vector (along events)
#' @return a matrix of Bernoulli likelihoods

mApplyBern <- function(X,prob) {
    t ( apply(X,1,function(x) dbern(x,prob)) )
}


#' Calculate categorical likelihood row-wise for categorical events
#' @param X a matrix of categorical events (row-wise)
#' @param prob the categorical probability matrix (rows along events, cols along event values)
#' @return a matrix of categorical likelihoods

mApplyCat <- function(X,prob) {
    t ( apply(X,1,function(x) dcat(x,prob)) )
}


#' Estimate log document probabilites given specific Bernoulli parameters
#' @param X a matrix of binary events (row-wise)
#' @param logprob the Bernoulli probability
#' 
#' @examples {
#'   X=matrix(c(0,1,0,0,0,0,1,0),2,4,byrow=TRUE) # two documents of length 4
#'   prob=c(.1,.2,.1,.1) # prob per index
#'   dput(mApplyBern(X,prob)) # likelihood for each index
#'   #structure(c(0.9, 0.9, 0.2, 0.8, 0.9, 0.1, 0.9, 0.9), .Dim = c(2L, 4L))
#'   dput(estLogPxBernGthetaJ(X,log(prob))) 
#'   # c(-1.92551945940758, -2.73644967562391)
#' }

estLogPxBernGthetaJ <- function (X, logprob) {
  # prob is a vector of size #words
  # X is a matrix (nrow=n, ncol=#words)
  apply(mApplyBern(X,exp(logprob)), 1, function(x) cumsum(log(x))[length(x)] )
}


#' Estimate log document probabilites given specific Categorical parameters
#' @param X a matrix of categorical events (row-wise)
#' @param logprob the Categorical probability
#' 
#' @examples {
#'   X=matrix(c(1,2,1,1,1,1,2,1),2,4,byrow=TRUE) # two documents of length 4
#'   prob=matrix(c(.9,.8,.9,.9,.1,.2,.1,.1),4,2) # prob per index
#'   dput(mApplyCat(X,prob)) # likelihood for each index
#'   #structure(c(0.9, 0.9, 0.2, 0.8, 0.9, 0.1, 0.9, 0.9), .Dim = c(2L, 4L))
#'   dput(estLogPxCatGthetaJ(X,log(prob))) 
#'   # c(-1.92551945940758, -2.73644967562391)
#' }

estLogPxCatGthetaJ <- function (X, logprob) {
  # prob is a matrix of size #words x #classes
  # X is a matrix (nrow=n, ncol=#words)
  apply(mApplyCat(X,exp(logprob)), 1, function(x) cumsum(log(x))[length(x)] )
}


#' Computes the cumulative sum in terms of logarithmic in- and output
#' Useful to avoid numerical underflow when summing products of probabilities
#' When using this function, one can sum sums of log probabilities
#' See also: http://goo.gl/aJopi
#' @param logx a vector of log numbers (need not be probabilities)
#' @return the log of the sum of the exponentiated input
#' @examples {
#'   x=c(1,2,3)
#'   exp(logsum(log(x)))
#'   # 6
#' }

logsum <- function(logx) {
  mypi=max(logx)
  mysum=0
  for (i in 1:length(logx)) {
    mysum = mysum + exp(logx[i]-mypi)
  }
  mypi + log(mysum)
}


#' objective function for mixture of binomials EM: 
#' @param DPS documents weighted by cluster priors
#' @return sum of log likelihood of documents

oFMixBinEM<- function (DPS) {
  sum(apply(DPS, 2, function(logx) logsum(logx)))
}
 

#' Assign final indices to data by maximum posterior value.
#' @param PjV1 Posterior matrix view 1 (by document).
#' @param PjV2 Posterior matrix view 2 (by document).
#' @return vector of indices for each data point.

assignIdxPerClMBinEM <- function(PjV1, PjV2) {
  apply( PjV1 + PjV2, 1, function(x) which.max(x) )
}


#' Agreement rate by maximum posterior values.
#' @param PjV1 Posterior matrix view 1 (by document).
#' @param PjV2 Posterior matrix view 2 (by document).
#' @return agreement rate.

agreementRateBinM  <- function(PjV1, PjV2) {
  v1 = apply( PjV1, 1, function(x) which.max(x) )
  v2 = apply( PjV2, 1, function(x) which.max(x) )
  sum(v1 == v2) / length(v1)
}
