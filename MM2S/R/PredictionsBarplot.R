# Deena M.A. Gendoo
# October 26, 2015
# Code to generate stacked barplot of MM2S predictions, per sample

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of this MM2S.

#################################################################################
#################################################################################

PredictionsBarplot<-function(InputMatrix,pdf_output,pdfheight,pdfwidth)
{
  colorscheme<-c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
  
  if(is.logical(pdf_output))
  {
    if(pdf_output==TRUE)
    {
      if((is.numeric(pdfheight))&&(is.numeric(pdfwidth)))
      {
      pdf(file="MM2S_Predictions_StackedBarplot.pdf",height=pdfheight,width=pdfwidth)
      par(c(2,1))
      barplot(t(InputMatrix),col=colorscheme,las=2,main="Human Subtype Predictions")
      legend("bottom",xpd = TRUE, inset=c(0,-0.1),fill = colorscheme,legend = colnames(InputMatrix),bty = "n",horiz = T)
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
  
  barplot(t(InputMatrix),col=colorscheme,las=2,main="Human Subtype Predictions")
  legend("bottom",xpd = TRUE, inset=c(0,-0.1),fill = colorscheme,legend = colnames(InputMatrix),horiz = TRUE)
}
