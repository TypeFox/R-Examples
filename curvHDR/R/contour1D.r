########## R function: contour1D ##########

# Obtains `contours' of a univariate function,
# for a given y-value and a grids of abscissa
# and ordinate values.

# Last changed: 14 APR 2009

contour1D <- function(yVal,xg,yg)
{  
   if ((yVal>max(yg))|(yVal<min(yg)))
      stop("yVal is outside the range of plausible values.")

   discrim <- diff(sign(c(yg-yVal)))
   discrim <- c(0,discrim)
   SolnInds <- (1:length(xg))[discrim!=0]
   xSolns <- (xg[SolnInds] + xg[SolnInds-1])/2

   numInts <- length(xSolns)/2
   lowLims <- xSolns[seq(1,(2*numInts-1),by=2)]
   uppLims <- xSolns[seq(2,(2*numInts),by=2)]

   return(lapply(1:length(lowLims),function(i)
           {return(c(lowLims[i],uppLims[i]))}))
}

########## End of contour1D ############
