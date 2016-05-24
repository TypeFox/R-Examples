########## R function: addIntervs ##########

# For adding a possibly disjoint set of
# intervals to an existing plot.

# Last changed: 17 JUL 2009

addIntervs <- function(intervList,stripWidth,polyCol="blue")
{
   # Convert list to a matrix:

   intervMat <- matrix(NA,length(intervList),2)
   for (i in 1:nrow(intervMat))
   {  
      intervMat[i,1] <- intervList[[i]][1]
      intervMat[i,2] <- intervList[[i]][2]
   }  
     
   # Ensure remaining intervals are properly sorted:

   orderVal <- order(intervMat[,1])
   if (length(orderVal)>1)
      intervMat <- intervMat[orderVal,]

   # Ensure that intervals are legals:

   if (!all(intervMat[,2]>=intervMat[,1]))
      stop("some intervals not legal.")

   # Ensure that intervals have correct monotonicity properties:

   if (!all(order(intervMat[,2])==1:nrow(intervMat)))
      stop("intervals are not monotonic.")

   # Set the matrix of remaining intervals.
   
   intsRemain <- intervMat
   indsToKeep <- 1:nrow(intsRemain)

   curPolLow <- 0 ; curPolUpp <- stripWidth
   while (!is.null(indsToKeep))
   {
      indsToKeep <- NULL
      numLeft <- nrow(intsRemain)

      curPoly <- list(x=c(intsRemain[1,],rev(intsRemain[1,])),
                      y=c(rep(curPolLow,2),rep(curPolUpp,2)))
      polygon(curPoly,col=polyCol)
      if (numLeft>1)
      {
         for (i in 2:numLeft)
         {
            if (intsRemain[i,1]>intsRemain[(i-1),2])
            {
               curPoly <- list(x=c(intsRemain[i,],rev(intsRemain[i,])),
                               y=c(rep(curPolLow,2),rep(curPolUpp,2)))
               polygon(curPoly,col=polyCol)
            }
            if (intsRemain[i,1]<=intsRemain[(i-1),2])
               indsToKeep <- c(indsToKeep,i)
         }
      }

      if (!is.null(indsToKeep))
         intsRemain <- intsRemain[indsToKeep,]

      if (length(indsToKeep)==1)
         intsRemain <- t(as.matrix(intsRemain))

      curPolLow <- curPolLow + 1.2*stripWidth
      curPolUpp <- curPolUpp + 1.2*stripWidth
   }
}

############ End of addIntervs ############

