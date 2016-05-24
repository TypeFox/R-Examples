# Version: 07-07-2013

# Changes:
#   Added the options "rejLine" and "alpha"

# Version: 25-05-2013, Daniel Fischer
# Changes:
#   Added the rej.lty parameter in order to change the line type of the rejection line
#   Added the rej.col parameter in order to change the color of the rejection line
#   Changed the default of crit to NULL, so that by default the lower part of the plot isn't plotted!
#   rej.para accepts now also vectors as entry and plots then several rejection lines.
#   Fixed a bug in the matrix input (couldn't handle data.frame input)

# 07-07-2013: Changed the Bonferroni slope (it depends now also on the alpha value
# 21-09-2013: Added the option for a legend
#             Added automatical colors, if not enough are given.

rejectionPlot <- function(X, lCol="red", xlim=NULL, crit=NULL, rejLine=NULL, alpha=0.01, rej.lty=c("dotted"), rej.col="black", incLegend=NULL){

  # Input check for the option 'crit'
    if(!is.null(crit)){
      crit <- match.arg(crit,c("distance","ratio"))
    } else {
      crit <- "NULL"
    }
  # Input check for the option 'rejLine'  
    if(!is.null(rejLine)) rejLine <- match.arg(rejLine,c("bh","bonferroni","simes"))
  # X-axis dimension setting
    if(is.null(xlim)) xlim <- c(0,1)
  # If the input is a data.frame bring it into matrix format
    if(is.data.frame(X)) X <- as.matrix(X)
  # Logical, was the input 'X' originally a matrix?  
    wasMatrix <- is.matrix(X) && dim(X)[1]>1
  # In case of vector input, transform it into a matrix with one row
    if(!is.matrix(X)) X <- t(as.matrix(X))
  # check if we have for each row in the input 'X' a different color
    if(length(lCol)!=dim(X)[1])
    {
      warning("Too less colors given, colors have automatically be assigned!")
      lCol <- 1:dim(X)[1]
    }

  # Get the required p-value grid (in order to calcualte the obs/exp plot - we pool the p-values over all given vectors, in order to calcualte the most dense one)
  temp <- unique(c(0,sort(as.vector(X)),1))
  sigTests <- matrix(NA,ncol=length(temp),nrow=dim(X)[1]+1)
  sigTests[1,] <- temp

  for(i in 1:length(temp))
  {
    temp <- (X <= sigTests[1,i])
    for(j in 2:dim(sigTests)[1])
    {
      sigTests[j,i] <- sum(temp[j-1,])
    }
  }
  for(i in 2:dim(sigTests)[1])
  {
    sigTests[i,] <- sigTests[i,]/dim(X)[2]
  }
  
  ratios <- matrix(NA,ncol=ncol(sigTests),nrow=nrow(sigTests)-1)
  distances <- matrix(NA,ncol=ncol(sigTests),nrow=nrow(sigTests)-1)
  for(i in 1:dim(ratios)[1])
  {
    ratios[i,] <- sigTests[i+1,]/sigTests[1,]
    distances[i,] <- sigTests[i+1,] - sigTests[1,]
  }
  ratios[1] <- 0
  
  ylim <- c(0,max(sigTests[,sum(sigTests[1,]<=xlim[2])]))
  
  if(crit=="NULL"){
    par(oma=c(2,0,1,0),
        mar=c(4,4,1,1))
  } else {
    nf <- layout(matrix(c(1,2),ncol=1),c(4,4), c(3,1), TRUE)
    par(oma=c(2,0,1,0),
        mar=c(4,4,1,0))
  }
  
  if (crit=="NULL"){
    plot(c(0,sigTests[1,]),c(0,sigTests[2,]),col=lCol[1],type="l",xlab="Expected Ratio",ylab="Observed Ratio",xlim=xlim,ylim=ylim)
  } else {
    plot(c(0,sigTests[1,]),c(0,sigTests[2,]),col=lCol[1],type="l",ylab="Observed Ratio",xlab=" ",xlim=xlim,ylim=ylim)
  }
  lines(c(-10,10),c(-10,10),type="l")
  if(wasMatrix)
  {
    for(i in 3:dim(sigTests)[1])
    {
      lines(c(0,sigTests[1,]),c(0,sigTests[i,]),type="l",col=lCol[i-1])
    }
  }

 # Add the different rejection lines, if requested:
  if(!is.null(rejLine))
  {
    if(rejLine=="bh"|rejLine=="simes")
    {
       for(i in 1:length(alpha)){
         lines(c(0,1),c(0,(1+alpha[i])/alpha[i]),lty=rej.lty, col=rej.col) 
       }
    } else if(rejLine=="bonferroni") {
       for(i in 1:length(alpha)){
         lines(c(0,1),c(0,dim(X)[2]/alpha[i]),lty=rej.lty, col=rej.col)
       }
    }
  }
 
# Add a possible Legend 
  if(!is.null(incLegend)){
    # legend(incLegend,legend=paste("Row",lCol), col=lCol, pch=20, lty=1,bg="white")#, trace=TRUE)
     legend(incLegend,legend=rownames(X), col=lCol, pch=20, lty=1,bg="white")#, trace=TRUE)
  }

 if(crit=="ratio")
 {
   plot(c(-1,2),c(1,1),type="l",xlim=xlim,ylim=c(min(ratios),max(ratios)),xlab="Expected Ratio",ylab="Ratio")
   for(i in 1:dim(ratios)[1])
   {
     lines(sigTests[1,],ratios[i,],col=lCol[i])
   }
  } else if(crit=="distance") {
     plot(c(-1,2),c(1,1),type="l",xlim=xlim,ylim=c(min(distances),max(distances)),xlab="Expected Ratio",ylab="Distance")
     for(i in 1:dim(ratios)[1])
     {
       lines(sigTests[1,],distances[i,],col=lCol[i])
     }
    lines(c(-10,10),c(0,0),lty="dotted")
  }
}