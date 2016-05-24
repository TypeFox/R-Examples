traceGeneProfiles <-
function(targetData, parentData, dataDescription=NULL,
			    targetColor = "grey", 
			    parentColor = "blue", onepage=TRUE)
{
  # In case of repeted measurements
  if(!is.null(dataDescription)){  
    
    TimePoints = unique(dataDescription)
    
    # Only one target gene to plotted
    if(is.null(nrow(targetData)))
    {
      MeanProfile = NULL
      for(i in 1:length(TimePoints))
      {
	MeanProfile = c(MeanProfile, 
		mean(targetData[dataDescription == TimePoints[i]]))
      }
      targetData = MeanProfile
      
    }else{
      MeanMatrix = matrix(0, nrow = nrow(targetData), ncol = length(TimePoints)) 

      for(i in 1:length(TimePoints))
      {

	for(j in 1:nrow(targetData))
	{
	  MeanValue = mean(targetData[j,dataDescription == TimePoints[i]])
	  MeanMatrix[j,i] = MeanValue
	}  
      }
      targetData = MeanMatrix
      
    }
  
    # Only one target gene to plotted
    if(is.null(nrow(parentData)))
    {
      MeanProfile = NULL
      for(i in 1:length(TimePoints))
      {
	MeanProfile = c(MeanProfile, 
		mean(parentData[dataDescription == TimePoints[i]]))
      }
      parentData = MeanProfile
      
    }else{
	MeanMatrix = matrix(0, nrow = nrow(parentData), ncol = length(TimePoints)) 

      for(i in 1:length(TimePoints))
      {

	for(j in 1:nrow(parentData))
	{
	  MeanValue = mean(parentData[j,dataDescription == TimePoints[i]])
	  MeanMatrix[j,i] = MeanValue
	}  
      }
      parentData = MeanMatrix
      
    }
 
  # End of if() dataDescription do exist
  }
  
  # Parent and target gene expression profiles will be represented here
  if(onepage)par(mfrow=c(2,1))

  # For limited the axes values
  MinY = min(rbind(targetData, parentData))
  MaxY = max(rbind(targetData, parentData))

  # Only one target gene is going to be plotted
  if(nrow(targetData)==1 || is.null(nrow(targetData)))
  {
      plot(1:length(targetData), targetData, type="l", ylim=c(MinY, MaxY),
      axes = FALSE,
      xlab = "Time point", ylab = "Expression value", 
      main = "Target gene expression profile", col = targetColor)
  }else{ # Several target genes are plotted

      plot(1:ncol(targetData), targetData[1,], type="l", ylim=c(MinY, MaxY),
      axes = FALSE,
      xlab = "Time point", ylab = "Expression value", 
      main = "Target gene expression profiles", col = targetColor)
      
      for (i in 2:nrow(targetData)){
	  lines(1:ncol(targetData), targetData[i,], col= targetColor, lwd = 1)
      }

  }

  axis(1)
  axis(2)
  
  #abline(v = CPpos, lty = "dashed")

  # Only one parentGene gene is going to be plotted
  if(nrow(parentData)==1 || is.null(nrow(parentData)))
  {
      plot(1:length(parentData), parentData, type="l", ylim=c(MinY, MaxY),
      axes = FALSE,
      xlab = "Time point", ylab = "Expression value", 
      main = "Parent gene expression profile", col = parentColor)
  }else{ # Several target genes are plotted

      plot(1:ncol(parentData), parentData[1,], type="l", ylim=c(MinY, MaxY),
      axes = FALSE,
      xlab = "Time point", ylab = "Expression value", 
      main = "Parent gene expression profiles", col = parentColor)
      
      for (i in 2:nrow(parentData)){
	  lines(1:ncol(parentData), parentData[i,], col= parentColor, lwd = 1)
      }

  }

  axis(1)
  axis(2)

  #abline(v = CPpos, lty = "dashed")

# End of function traceGeneProfiles
}
