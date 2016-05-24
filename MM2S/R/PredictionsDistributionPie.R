# Deena M.A. Gendoo
# October 26, 2015
# Code to generate stacked barplot of MM2S predictions, per sample

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of this MM2S.

#################################################################################
#################################################################################

PredictionsDistributionPie<-function(InputMatrix,pdf_output,pdfheight,pdfwidth)
{
  colorscheme<-c("#fccde5","#ccebc5","#80b1d3","#bc80bd","#fb8072")
  
  if(is.logical(pdf_output))
  {
    if(pdf_output==TRUE)
    {
      if((is.numeric(pdfheight))&&(is.numeric(pdfwidth)))
      {
      pdf(file="MM2S_OverallPredictions_PieDist.pdf",height=pdfheight,width=pdfwidth)
      pie(c(table(InputMatrix$MM2S_Subtype[,2])),main="Overall Predictions",labelcex = 0.6,col=colorscheme)
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
  
  pie(c(table(InputMatrix$MM2S_Subtype[,2])),main="Overall Predictions",labelcex = 0.6,col=colorscheme)
}
