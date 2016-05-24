plotSegmentModel.postDist <-
function(edgesPostDist, CPpos, parentNames=NULL, targetName=NULL, edgesThreshold=0.5, onepage=TRUE){
  nbSegs=length(CPpos)-1
  if(onepage)par(mfrow=c(1,nbSegs))
  if(is.null(parentNames)) parentNames=c(1:(dim(as.matrix(edgesPostDist))[2]))
  
  for(i in 1:nbSegs){
      
      SelectedParents = length(which(edgesPostDist[i,] > edgesThreshold ))
      if(max(nchar(parentNames))>2){
        graphic.las=2
        par(mar=c(6,4,4,2),mgp=c(3, 0.6, 0))
      }else{ graphic.las=1}
      
      barplot(edgesPostDist[i,],names.arg=substr(parentNames,1,6),
	      #xlab=paste("Parent genes"),
              ylab="Estimated posterior probability",
	      main=paste("Regulatory model for target gene:", targetName, "\n", "Temporal segment #", i, ":  [",CPpos[i],",",CPpos[i+1]-1 ,"] \n"),
      	      lwd=2,col="blue",ylim=c(0,1), las= graphic.las,font.lab=2) 
      if(max(nchar(parentNames))>2)par(mgp=c(4.7,0.1,0))
      title(xlab=paste("Parent genes \n ","# of selected edge(s):", SelectedParents),font.lab=2)
      par(mgp=c(3, 1, 0))
  
      abline(h = edgesThreshold, lty = "dashed", col = "grey")
    
    }
 }
