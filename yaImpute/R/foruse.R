# Takes a "yai" object and returns a "for-use" table.
#
# When method=="kth"
# By default (kth=NULL), only the best pick for an observation is listed in the
# use column. For a reference observations, it is always selected to represent
# itself. However, when kth is not NULL, the kth neighbor is reported.
#
# When method=="random" or "randomWeighted"
# For each target, kth is selected at random from the k neighbors. When
# "randomWeighted" is used, 1/(1+d) is used as a probability weight factor
# in selecting the kth neighbor, where d is the distance.
# In this case kth is normally NULL, and set to k for the object. When kth is
# not null, it is used as the upper limit on the number of neighbors to consider. 
#
# When targetsOnly is true, reporting of references is not done.
#
# The rowname is a target id, first column is the reference id that
# is select to represent the target, the second is the corresponding
# distance. Every reference is included as a target as well. In those
# cases, the rowname and the use value are the same and the distance
# is zero.

foruse = function (object,kth=NULL,method="kth",targetsOnly=FALSE)
{
   if (class(object) != "yai") stop ("class must be yai")
   valid=c("kth","random","randomWeighted")
   if (is.na(match(method,valid))) stop (paste("method must be one of",paste(valid,collapse=", ")))
   
   if (method != "kth" && is.null(kth)) kth=object$k
   
   if (!is.null(kth))
   {
      if (kth>object$k) kth=object$k
      if (kth<1)        kth=NULL
   }
   if (!is.null(object$neiIdsRefs) && !targetsOnly)
   {
      if (is.null(kth))
         fu1=data.frame(use=rownames(object$neiIdsRefs),
                        dist=rep(0,nrow(object$neiIdsRefs)),
                        stringsAsFactors = FALSE)
      else
      {
         if (method != "kth")
         {
           ans = lapply(1:nrow(object$neiIdsRefs),
              function (x,kth,Ids,Dst,method) 
              {
                k = sample.int(kth,1,prob=if (method=="random") NULL else 1/(Dst[x,]+1))
                list(Ids[x,k], Dst[x,k])
              }, kth, object$neiIdsRefs, object$neiDstRefs, method)
           fu1=data.frame(use=unlist(lapply(ans,function (x) x[[1]])),
                          dist=unlist(lapply(ans,function (x) x[[2]])),
                          stringsAsFactors = FALSE)
         }
         else
         {       
           if (is.null(kth)) kth=1
           fu1=data.frame(use=object$neiIdsRefs[,kth],
                          dist=object$neiDstRefs[,kth],
                          stringsAsFactors = FALSE)
         }
      }
      rownames(fu1)=rownames(object$neiIdsRefs)
   }

   else fu1=NULL

   if (!is.null(object$neiIdsTrgs))
   {
      if (method != "kth")
      {
        ans = lapply(1:nrow(object$neiIdsTrgs),
           function (x,kth,Ids,Dst,method) 
           {
             k = sample(kth,1,prob=if (method=="random") NULL else 1/(Dst[x,]+1))
             list(Ids[x,k], Dst[x,k])
           }, kth, object$neiIdsTrgs, object$neiDstTrgs, method)
        fu2=data.frame(use=unlist(lapply(ans,function (x) x[[1]])),
                       dist=unlist(lapply(ans,function (x) x[[2]])),
                       stringsAsFactors = FALSE)
      }
      else
      { 
        if (is.null(kth)) kth=1
        fu2=data.frame(use=object$neiIdsTrgs[,kth],
                       dist=object$neiDstTrgs[,kth],
                       stringsAsFactors = FALSE)
      }                      
      rownames(fu2)=rownames(object$neiIdsTrgs)
   }
   else fu2=NULL
   if (is.null(fu1) & is.null(fu2)) return (NULL)
   if      (is.null(fu1)) ans = fu2
   else if (is.null(fu2)) ans = fu1
   else                   ans = rbind(fu1,fu2)
   class(ans)=c("data.frame","foruse.yaImpute")
   ans
}
