TopList <- function(object, topnum=10, sorted.by = c("RValue","PostMean","MLE","PVal"))  {
  ####  Example usage:
  ####  toplist(rvs,topnum=15,sorted.by="MLE")

  sorted.by <- match.arg(sorted.by)
  om <- object$main
  switch(sorted.by,
         RValue = head(om,topnum),
         PostMean = om[order(om$PM.rank)[1:topnum],],
         MLE = om[order(om$MLE.rank)[1:topnum],],
         PVal = om[order(om$PVal.rank)[1:topnum],])
}
