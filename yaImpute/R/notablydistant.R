# Provides a list of notably distant targets.
#
# Arguments:
#   object is of class yai (built using function yai), it
#          must contain reference distances if threshold is null.
#   threshold is a threshold value used in the calculations, if
#             NULL it is computed using the p value.
#   p is the percentile point in the log normal distribution used
#     to compute the threshold when it is null.
#   method="quantile": the threshold is computed using function "quantile" 
#     to pick the (1-p)th percentile point in the set of distances.
#     "distribution": the threshold is based on a percentile points
#     of a assumed distribution.
#     and otherwise the "distribution" assumptions are followed.
# Value:
#  List of two data frames that contain 1) the references that are notably
#  distant from other references, and 2) the targets that are notably distant
#  from the references, the threshold used and lastly the method used.

notablyDistant  <-  function (object,kth=1,threshold=NULL,p=0.01,
                              method="distribution")
{
   if (missing(object)) stop ("object required.")
   if (class(object) != "yai") stop ("class must be yai")
   if (kth>object$k) kth <- object$k
   if (kth<1)        kth <- 1
   if (is.null(threshold))
   {
      threshold <- NA
      if (is.null(object$neiDstRefs)) stop ("distances among references are required when threshold is NULL")
      
      if (method=="distribution")
      {      
        # use the beta disrtibution, distances are 0<=d<=1
        if (object$method %in% c("randomForest","random")) 
        {
           m <- mean(object$neiDstRefs[,kth])
           ss <- var(object$neiDstRefs[,kth])
           if (!is.nan(ss) & !is.nan(m))
           {
              v <- m*((m*(1-m)/ss)-1)
              w <- (1-m)*((m*(1-m)/ss)-1)
              threshold <- qbeta(p,v,w,lower.tail=FALSE)
           }
        }
        else # use the lognormal distribution, distances are 0<=d
        {
           zeros <- object$neiDstRefs[,kth]<=0
           if (sum(zeros)==0) obs <- log(object$neiDstRefs[,kth])
           else
           {
              smz <- min(object$neiDstRefs[!zeros,kth])
              obs <- object$neiDstRefs[,kth]
              obs[zeros] <- smz*.5
              obs <- log(obs)
              warning ("when computing threshold, ",sum(zeros)," zero distances of ",
                        length(obs)," references were set to ",format(smz*.5))
           }
           m <- mean(obs)
           s <- sd(obs)
           if (!is.nan(s) & !is.nan(m)) threshold <- exp(s*qnorm(p, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)+m)
        }
        if (is.nan(threshold))
        {
           threshold <- Inf
           warning ("threshold can not be computed, set to Inf")
        }
      }
      else
      {
        if (method != "quantile") 
        {
          method="quantile"
          warning("method set to quantile")
        }
        threshold <- quantile(object$neiDstRefs[,kth],probs=1-p)
      } 
   }
   findNDist <- function (ids,dst,names,threshold)
   {
      out <- data.frame(use=ids,dist=dst,row.names=names)
      out <- out[out[,2]>threshold,]
      if (nrow(out)>1)
      {
        ix  <- sort(out[,2],decreasing=TRUE, index.return=TRUE)$ix
        out <- out[ix,]
      }
      out
   }

   if (!is.null(object$neiDstRefs))
      distRefs <- findNDist(object$neiIdsRefs[,kth],object$neiDstRefs[,kth],rownames(object$neiIdsRefs),threshold)
   else distRefs=NULL

   if (!is.null(object$neiIdsTrgs))
      distTrgs <- findNDist(object$neiIdsTrgs[,kth],object$neiDstTrgs[,kth],rownames(object$neiIdsTrgs),threshold)
   else distTrgs=NULL

   list(notablyDistantRefs=distRefs, notablyDistantTrgs=distTrgs, threshold=threshold, method=method)
}
