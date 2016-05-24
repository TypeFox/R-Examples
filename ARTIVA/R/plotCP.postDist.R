plotCP.postDist <-
function(CPpostDist,targetName=NULL, onepage=TRUE, color1="green", color2="black", estimatedCPpos=NULL )
{
  CPnumberPostDist=CPpostDist$CPnumber
  CPpositionPostDist=CPpostDist$CPposition
  if(onepage)par(mfrow=c(1,2))
  
	MaxValue = which.max(CPnumberPostDist) - 1
        barplot(CPnumberPostDist,names.arg=0:(length(CPnumberPostDist)-1),
	    xlab= paste("Changepoint number\n", "( max. value =", MaxValue, ")"),
	    ylab="Estimated posterior probability",
	    main=paste("Number of changepoint\n", "- target gene:", targetName, "-"),
	    lwd=2,col=color1,ylim=c(0,1))
    abline(v = MaxValue+1, lty = "dashed", col = "red")
    
    plot(1:length(CPpositionPostDist),CPpositionPostDist,type="h",
    xlab=paste("Time point\n", "( # of selected changepoint(s) =", MaxValue, ")"),
    ylab="Estimated posterior probability",
    main=paste("Changepoint position\n", "- target gene:", targetName, "-"),
    lwd=2,col=color2,ylim=c(0,1))
  
    if(!is.null(estimatedCPpos)){
      VecMaxValue = estimatedCPpos[-c(1,length(estimatedCPpos))]
      abline(v = VecMaxValue, lty = "dashed", col = "red")
    }else{
      if(MaxValue > 0)
        {
          VecMaxValue = order(CPpositionPostDist, decreasing = TRUE)[1:MaxValue]
          abline(v = VecMaxValue, lty = "dashed", col = "red")
        }
    }
  
}
