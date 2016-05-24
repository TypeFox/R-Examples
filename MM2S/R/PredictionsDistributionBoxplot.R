# Deena M.A. Gendoo
# October 26, 2015
# Code to generate stacked barplot of MM2S predictions, per sample

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of this MM2S.

#################################################################################
#################################################################################

PredictionsDistributionBoxplot<-function(InputMatrix,pdf_output,pdfheight,pdfwidth)
{
  colorscheme<-c("#fccde5","#ccebc5","#80b1d3","#ffffb3","#fb8072")
  
  if(is.logical(pdf_output))
  {
    
    TrueCounts<-apply(InputMatrix$Predictions,2,function(x){sum(as.numeric(x)>0)})
    if(pdf_output==TRUE)
    {
      if((is.numeric(pdfheight))&&(is.numeric(pdfwidth)))
      {
        
      
        
      pdf(file="MM2S_OverallPredictions_BoxplotDist.pdf",height=pdfheight,width=pdfwidth)
        boxplot(InputMatrix$Predictions,las=2,ylab="Prediction Strength (%)",
                names=paste(colnames(InputMatrix$Predictions), " (n=", TrueCounts,")",sep=""),col=colorscheme,cex.axis=0.5,
                main=paste("Distribution of MB Subtype Calls Across ",nrow(InputMatrix$Predictions)," Samples",sep=""))
      dev.off()
      }
      else
      {message("PDF dimensions must be numeric")
      stop()}
    }
  }
  else {
    message("TRUE or FALSE needed for PDF output")
    stop()}
  
  boxplot(InputMatrix$Predictions,las=2,ylab="Prediction Strength (%)",
          names=paste(colnames(InputMatrix$Predictions), " (n=", TrueCounts,")",sep=""),col=colorscheme,cex.axis=0.5,
          main=paste("Distribution of MB Subtype Calls Across ",nrow(InputMatrix$Predictions)," Samples",sep=""))
}
