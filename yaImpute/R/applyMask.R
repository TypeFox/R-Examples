applyMask <- function
(object, refGroups=NULL, trgGroups=NULL, method="removeWhenCommon", k=1)
{
  if (class(object) != "yai") stop("object must be of class yai")
  valid <- c("removeWhenCommon","keepWhenCommon")
  if (is.na(match(method,valid))) stop (paste("method must be one of",paste(valid,collapse<-", ")))
  if (is.null(refGroups) | is.null(trgGroups)) stop("refGroups and trgGroups must be defined")
  if (k >= object$k) stop("new value of k (",k,") must be less than old value (",object$k,")")

  object$call <- match.call()
  refGrp <- refGroups[match(object$neiIdsTrgs,rownames(object$xRefs))]
  lrefGrp <- if (method == "removeWhenCommon") refGrp != trgGroups else refGrp == trgGroups
  dim(lrefGrp) <- dim(object$neiIdsTrgs) 

  # tvec is an offset in the storage of neiIdsTrgs and neiDstTrgs. At this point
  # The kth member is the offset of the first row of the kth column
  tvec <- 0:(ncol(lrefGrp)-1) * nrow(lrefGrp)
  
  # ans is the value of tvec corresponding the the columns to keep for each row.
  ans <- apply(lrefGrp,1,function(x,tvec,k) tvec[x][1:k],tvec,k)
  
  # if k>1, we need to reorganize ans and delete the dimensions so it is a vector.
  if (k>1) 
  {
    ans <- t(ans) 
    dim(ans) <- NULL
  }

  # now add the row numbers to ans...to get the final offsets. 
  ans <- rep(1:nrow(lrefGrp),k) + ans
  
  rnB <- rownames(object$neiIdsTrgs)
  cnI <- colnames(object$neiIdsTrgs)[1:k]
  cnD <- colnames(object$neiDstTrgs)[1:k]
  object$neiIdsTrgs <- object$neiIdsTrg[ans]
  object$neiDstTrgs <- object$neiDstTrg[ans]
  dim (object$neiIdsTrgs) <- c(nrow(lrefGrp),k)
  dim (object$neiDstTrgs) <- c(nrow(lrefGrp),k)
  rownames(object$neiIdsTrgs) <- rnB
  rownames(object$neiDstTrgs) <- rnB 
  colnames(object$neiIdsTrgs) <- cnI
  colnames(object$neiDstTrgs) <- cnD

  object$k <- k
  object$noRefs <- TRUE
  object$neiIdsRefs <- NULL
  object$neiDstRefs <- NULL
 
  object
}

