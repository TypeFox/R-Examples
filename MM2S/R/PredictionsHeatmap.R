# Deena M.A. Gendoo
# October 29, 2014
# Code to generate heatmap of MM2S predictions

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of this MM2S.

#################################################################################
#################################################################################

PredictionsHeatmap<-function(InputMatrix,pdf_output,pdfheight,pdfwidth)
{
  colorscheme<-c("#ffffff","#ffffcc","#c7e9b4","#41b6c4","#2c7fb8","#253494") 
  
  if(is.logical(pdf_output))
  {
    if(pdf_output==TRUE)
    {
      if((is.numeric(pdfheight))&&(is.numeric(pdfwidth)))
      {
      pheatmap(InputMatrix,color=colorscheme,cluster_rows=FALSE, cluster_cols=FALSE,cexRow=0.5,cexCol=0.5,border_color="black",
               scale="none",cellwidth=25,cellheight=15,legend=T,legend_breaks=c(0,20,40,60,80,100),
               legend_labels=c("0", "20", "40", "60", "80", "100"),main="Human Subtype Predictions",
               filename="MM2S_Predictions_Heatmap.pdf",height=pdfheight,width=pdfwidth)
      }
      else
      {message("PDF dimensions must be numeric")
      stop()}
    }
  }
  else {
    message("TRUE or FALSE needed for PDF output")
    stop()}
  
  pheatmap(InputMatrix,color=colorscheme,cluster_rows=FALSE, cluster_cols=FALSE,cexRow=0.5,cexCol=0.5,border_color="black",
           scale="none",cellwidth=25,cellheight=15,legend=T,main="Subtype Predictions")

}
