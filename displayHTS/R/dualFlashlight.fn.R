dualFlashlight.fn <-
function(result.df, wellName="WELL_USAGE", x.name="mean", y.name="ssmd",
         sampleName="Sample", sampleColor="black", 
         controls = NULL, controlColors = NULL, 
         xlab="Average Fold Change", ylab="SSMD",
         main="Dual-Flashlight Plot", x.legend=NA, y.legend=NA, 
         cex.point=1, cex.legend = 0.8, 
         xat=NULL, xMark=NULL, yat=NULL, yMark=NULL, xLines=NULL, yLines=NULL )
{
#*****************************************************************************#
# author: Xiaohua Douglas Zhang & Zhaozhi Zhang, 2012                         #
# function to draw dual-flashlight plot, volcano plot &plate correlation plot #
# Input:                                                                      #
#   result.df: results or data including at least three columns for well name,#
#         x (e.g., average fold change in log2) and y (e.g.,SSMD), respectively#
#   wellName: data column name indicating well names                          #
#   x.name: column name indicating values for x-axis                          #
#   y.name: column name indicating values for y-axis                          #
#   sampleName: name for the well type indicating sample wells                #
#   sampleColors: color for sample wells                                      #
#   controls: a vector including controls to be shown in the plot             #
#   controlColors: a vector including the color of the controls to be shown   #
#         in the plot. It must have the same length as 'controls'             #
#   xlab, ylab, main: same as those in the internal function 'plot'           #
#   x.legend, y.legend: the x and y co-ordinates to be used to position the   #
#         legend.                                                             #
#   cex.point: define the size of points in the plot                          #
#   cex.legend: define the size of legend in the plot                         #
#   xat: the position of x-axis to be labeled                                 #
#   xMark: the labels of x-axis corresponding to 'xat'                        #
#   xLines: x-values indicating positions of vertical grey lines to be drawn  #
#   yLines: y-values indicating positions of horizontal grey lines to be drawn#
# References:
#   Zhang XHD, Zhang ZZ. 2012. displayHTS: a R package displaying data and
#        results from high-throughput screening experiments (submitted).
#   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical Experimental
#        Design and Data Analysis for Genome-scale RNAi Research.
#        Cambridge University Press, Cambridge, UK.
#   Zhang XHD. 2009. A method effectively comparing gene effects in multiple
#        conditions in RNAi and expression profiling research.
#        Pharmacogenomics 10(3):345-358.
#   Zhang XHD. 2010. Assessing the size of gene or RNAi effects in multi-factor
#        high-throughput experiments. Pharmacogenomics 11(2): 199 - 213.   
#   Zhang XHD, Santini F, Lacson R, Marine SD, Wu Q, Benetti L, Yang R,
#        McCampbell A, Berger JP, Toolan DM, Stec EM, Holder DJ, Soper KA,
#        Heyse JF and Ferrer M. 2011. cSSMD: Assessing collective activity of
#        multiple siRNAs in genome-scale RNAi screens.
#        Bioinformatics 27(20): 2775-2781.  
# See Also:
#   imageDesign.fn, imageIntensity.fn, plateWellSeries.fn
# Example:                                                                    #
#  # for dual-flashlight plot
#  data("HTSresults", package = "displayHTS")
#  par( mfrow=c(1, 1) )
#  dualFlashlight.fn(HTSresults, wellName="WELL_USAGE", x.name="mean",
#                    y.name="ssmd", sampleName="Sample", sampleColor="black", 
#                    controls = c("negCTRL", "posCTRL1", "mock1"),
#                    controlColors = c("green", "red", "lightblue"), 
#                    xlab="Average Fold Change", ylab="SSMD",
#                    main="Dual-Flashlight Plot", x.legend=NA, y.legend=NA, 
#                    cex.point=1, cex.legend = 0.8,
#                    xat=log2( c(1/8, 1/4, 1/2, 1, 2, 4, 8) ), 
#                    xMark=c("1/8", "1/4", "1/2", "1", "2", "4", "8"),
#                    xLines=log2(c(1/4, 1/2 ,1, 2, 4)),
#                    yLines=c(-5, -3, -2, -1, 0, 1, 2, 3, 5 ) )
#  # for volcano plot
#  result.df=cbind(HTSresults,"neg.log10.pval"=-log10(HTSresults[,"p.value"]))
#  dualFlashlight.fn(result.df, wellName="WELL_USAGE", x.name="mean",
#                    y.name="neg.log10.pval",
#                    sampleName="Sample", sampleColor="black", 
#                    controls = c("negCTRL", "posCTRL1", "mock1"),
#                    controlColors = c("green", "red", "lightblue"), 
#                    xlab="Average Fold Change", ylab="p-value",
#                    main="Volcano Plot", x.legend=NA, y.legend=-log10(0.06), 
#                    cex.point=1, cex.legend = 0.8,
#                    xat=log2( c(1/8, 1/4, 1/2, 1, 2, 4, 8) ), 
#                    xMark=c("1/8", "1/4", "1/2", "1", "2", "4", "8"),
#                    yat=-log10( c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1) ), 
#                    yMark=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
#                    xLines=log2(c(1/4, 1/2 ,1, 2, 4)),
#                    yLines=-log10( c( 0.001, 0.01, 0.05) ) )
#  # plate pair correlation plot
#  data("HTSdataSort", package = "displayHTS")  
#  data.df= cbind(HTSdataSort[1:384,], HTSdataSort[384+1:384,])
#  names(data.df)=
#    c("SOBARCODE.1", "BARCODE.1", "XPOS.1", "YPOS.1", "WELL_USAGE.1",   
#      "Compound.1", "Intensity.1", "log2Intensity.1",
#      "SOBARCODE.2", "BARCODE.2", "XPOS.2", "YPOS.2", "WELL_USAGE.2",   
#      "Compound.2", "Intensity.2", "log2Intensity.2")
#  dualFlashlight.fn(data.df, wellName="WELL_USAGE.1", x.name="log2Intensity.1",
#                    y.name="log2Intensity.2", 
#                    sampleName="Sample", sampleColor="black", 
#                    controls = c("negCTRL", "posCTRL1", "mock1"),
#                    controlColors = c("green", "red", "lightblue"), 
#                    xlab="log2 intensity in plate 1",
#                    ylab="log2 intensity in plate 2",
#                    main="Plate Pair Correlation Plot", x.legend=NA,
#                    y.legend=NA, cex.point=1, cex.legend = 0.8 )
#  abline(0,1)
#*****************************************************************************#

  condt = result.df[,wellName] == sampleName
  Legends = sampleName
  if( !is.null(controls) ) {
    Legends = c(sampleName, controls)
    col.legend = c(sampleColor, controlColors)
    for( i in 1:length(controls) ) {
      condt = condt | result.df[,wellName] == controls[i]
    }
  }
  
  if( is.na(x.legend) ) x.legend = min(result.df[condt, x.name], na.rm=T)
  if( is.na(y.legend) ) y.legend = max(result.df[condt, y.name], na.rm=T)

  xRange = range(result.df[condt, x.name], na.rm=T)
  yRange = range(result.df[condt, y.name], na.rm=T)
  plot( xRange, yRange, type="n", xlab=xlab, ylab=ylab,
       cex=cex.point, axes=F, main=main )
  if( is.null(xat) ) axis(1) else axis(1, at=xat, labels=xMark )
  if( is.null(yat) ) axis(2, las=2) else axis(2, at=yat, labels=yMark, las=2 )
  box()
  condt = result.df[, wellName ] == sampleName
  points(result.df[condt, x.name],
         result.df[condt, y.name], col=sampleColor, cex=cex.point)

  if( !is.null(controls) ) {
    for( i in 1:length(controls) ) {
     condt = result.df[, wellName ] == controls[i]
     points(result.df[condt, x.name],
            result.df[condt, y.name], col=controlColors[i], cex=cex.point)
    }
  }

  if( !is.null(xLines) )
    for(a in xLines) lines(rep(a, 2), c(-10000, 10000), col="grey")
  if( !is.null(yLines) )
    for(a in yLines) lines(c(-10000, 10000), rep(a, 2), col="grey")

  legend(x.legend, y.legend, legend=Legends, col=col.legend, cex=cex.legend,
         pch=1, bg="white")
}
