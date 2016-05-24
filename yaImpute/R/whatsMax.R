# For each row, finds the column that has the maximum value. Returns
# a data frame with two columns: The first is the column name corresponding
# to the column of maximum value and the second is the correspond maximum.
# The first column is convereted to a factor.
#
# if the max is zero, the maxCol is identified as "zero".
# if there are too many "factors" (too many meaning over "nbig") in column 1,
# the lowest maxs are all termed "other".
#
# Intended use is to transform community ecology data for use in yai
# where method is randomForest.

whatsMax <- function (x,nbig=30)
{
   num <- vector("logical",ncol(x))
   for (i in 1:ncol(x)) num[i] <- is.numeric(x[,i])
   if (sum(as.numeric(num))==0) stop ("no numeric columns")
   ag <- deparse(substitute(x))
   orgRows <- nrow(x)
   x <- na.omit(x)
   if (nrow(x) != orgRows) warning (orgRows-nrow(x)," rows had missing data and were deleted")
   n <- colnames(x)
   ans <- apply(x,1,function (x1) which.max(x1))
   nm  <- lapply(ans,function (i,n) n[i], colnames(x))
   ans <- apply(x,1,function (x1) x1[which.max(x1)])
   ans <- as.data.frame(list(a=unlist(nm),unlist(ans)),stringsAsFactors=FALSE)
   rownames(ans) <- rownames(x)
   colnames(ans) <- paste(ag,c("maxCol","maxVal"),sep=".")
   ans[ans[,2]==0,1] <- "zero"
   nf <- length(levels(as.factor(ans[,1])))
   if (nf > nbig)
   {
      tops <- names(sort(tapply(ans[,2],as.factor(ans[,1]),sum),decreasing = TRUE))
      ans[ans[,1] %in% tops[nbig:length(tops)],1] <- "other"
   }
   ans[,1] <- as.factor(ans[,1])
   ans[,2] <- as.numeric(ans[,2])
   ans
}
