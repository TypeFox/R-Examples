plateWellSeries.fn <-
function(data.df, intensityName=NA, plateName="BARCODE", wellName="WELL_USAGE",
         rowName, colName, show.wellTypes=NULL, order.wellTypes=NULL,
         color.wells=NULL, pch.wells=NULL, ppf=12, byRow=T, yRange=NULL,
         cex.point=0.25, cex.legend=0.3, x.legend=NA, y.legend=NA,
         main=NA, xlab=NA, ylab=NA)
{
#*****************************************************************************#
# author: Xiaohua Douglas Zhang & Zhaozhi Zhang, 2012                         #
# function to display the image of well designs in a plate                    #
# Input:                                                                      #
#   data.df: data for all plates including at least five columns for          #
#         intensity, plate names, well names, rows and columns, respectively  #
#   wellName: data column name indicating well names                          #
#   rowName:  data column name indicating row numbers in a plate              #
#   colName:  data column name indicating column numbers in a plate           #
#   show.wellTypes: a vector of well types to be shown in the well-series plot#
#   order.wellTypes: a vector of numbers to indicate the order of well types  #
#         corresponding to 'show.wellTypes' in the well-series plot           # 
#   color.wells:  a vector indicating the colors of well types corresponding  #
#         to 'show.wellTypes' in the well-series plot                         #
#   pch.wells:  a vector indicating the point types of well types             #
#         corresponding to 'show.wellTypes' in the well-series plot           #
#   ppf: number of plates per figure to be shown in the well-series plot      #
#   byRow: indicating whether the wells in a plate should be shown by row     #
#         or column                                                           #
#   yRange: define the range of y-axis in the well-series plot                #
#   cex.point: define the size of points in the well-series plot              #
#   cex.legend: define the size of legend in the well-series plot             #
# References:
#   Zhang XHD, Zhang ZZ. 2012. displayHTS: a R package displaying data and
#        results from high-throughput screening experiments (submitted).  
#   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical Experimental
#        Design and Data Analysis for Genome-scale RNAi Research.
#        Cambridge University Press, Cambridge, UK.
#   Zhang XHD, Yang XC, Chung N, Gates AT, Stec EM, Kunapuli P, Holder DJ,
#        Ferrer M, Espeseth AS. 2006. Robust statistical methods for hit
#        selection in RNA interference high throughput screening experiments.
#        Pharmacogenomics 7 (3) 299-309.
# See Also:
#   image.design.fn, image.intensity.fn, dualFlashlight.fn
# Example:                                                                    #
#   data("HTSdataSort", package = "displayHTS")
#   wells = as.character(unique(HTSdataSort[, "WELL_USAGE"])); wells
#   colors = c("black",  "yellow", "grey", "blue", "skyblue", "green", "red")
#   orders=c(1, 3, 2, 4, 5, 7, 6)
#   # by row
#   par( mfrow=c(2,1) )                                                        
#   plateWellSeries.fn(data.df = HTSdataSort, intensityName="log2Intensity",
#                      plateName="BARCODE", wellName="WELL_USAGE",              
#                      rowName="XPOS", colName="YPOS", show.wellTypes=wells,
#                      order.wellTypes=orders, color.wells=colors,
#                      pch.wells=rep(1, 7), ppf=6, byRow=T,  
#                      yRange=NULL, cex.point=0.25,cex.legend=0.3) 
#   # by column
#   par( mfrow=c(2,1) )                                                        
#   plateWellSeries.fn(data.df = HTSdataSort, intensityName="log2Intensity",
#                      plateName="BARCODE", wellName="WELL_USAGE",              
#                      rowName="XPOS", colName="YPOS", show.wellTypes=wells,
#                      order.wellTypes=orders, color.wells=colors,
#                      pch.wells=rep(1, 7), ppf=6, byRow=F,  
#                      yRange=NULL, cex.point=0.25,cex.legend=0.3) 
#   # display hits
#   data("HTSresults", package = "displayHTS")
#   condtSample = HTSresults[, "WELL_USAGE"] == "Sample"
#   condtUp = HTSresults[,"ssmd"] >= 1 & HTSresults[,"mean"] >= log2(1.2)
#   condtDown = HTSresults[,"ssmd"] <= -1 & HTSresults[,"mean"] <= -log2(1.2)
#   sum(condtSample & (condtUp | condtDown) )/sum(condtSample)
#   hit.vec = as.character(HTSresults[, "WELL_USAGE"])
#   hit.vec[ condtSample & condtUp ] = "up-hit"
#   hit.vec[ condtSample & condtDown ] = "down-hit"
#   hit.vec[ condtSample & !condtUp & !condtDown] = "non-hit"
#   result.df = cbind(HTSresults, "hitResult"=hit.vec)
#   wells = as.character(unique(result.df[, "hitResult"])); wells
#   orders = c(1, 3, 4, 6, 7, 8, 9, 2, 5) 
#   colors = c("black",  "green", "yellow", "red",
#              "grey", "purple1", "purple2", "lightblue", "purple3")
#   par(mfrow=c(1,1))
#   plateWellSeries.fn(data.df = result.df, intensityName="mean",       
#                      plateName="SOBARCODE", wellName="hitResult", 
#                      rowName="XPOS", colName="YPOS", show.wellTypes=wells,
#                      order.wellTypes=orders, color.wells=colors,
#                      pch.wells=rep(1, 7), ppf=6, byRow=F,  
#                      yRange=NULL, cex.point=0.5,cex.legend=0.55,
#                      y.legend=-0.5) 
#*****************************************************************************#

  plateWelltoX.fn = function(data.df, nRow=16, nCol=24, byRow=T) {
  #*************************************************************************
  # function to generate value in x-axis based on plate number and
  #   well positions in a plate.
  # data.df: consists of three columns: plate names or index, rows and columns
  # Author: Xiaohua Douglas Zhang, 2005
  #*************************************************************************
    plate.vec = data.df[,1]
    platesUniq = unique(plate.vec)
    plateOrder.vec = rep( NA, length(plate.vec) )
    
    for( i in 1:length(platesUniq) )
      plateOrder.vec[plate.vec==platesUniq[i]] = i
    if( byRow ) {
      x.vec = (plateOrder.vec-1)*nRow*nCol + (data.df[,2]-1)*nCol + data.df[,3]
    } else {
      x.vec = (plateOrder.vec-1)*nRow*nCol + (data.df[,3]-1)*nRow + data.df[,2]
    }    
    data.frame(x=x.vec, plateOrder = plateOrder.vec)
  }

  inten.vec0 = data.df[, intensityName]
  plate.vec = as.vector(data.df[, plateName])
  plates = unique( plate.vec )
  nPlate = length(plates)

  nRow=max(data.df[,rowName]); nCol=max(data.df[,colName]); nWell=nRow*nCol
 
  x.df=cbind(plateWelltoX.fn(data.df[,c(plateName, rowName, colName)],
                             nRow=nRow, nCol=nCol, byRow=byRow),
             "wellType"=data.df[, wellName], "plateName"= data.df[,plateName])
  for(k in 1:ceiling(nPlate/ppf) ) { 
    wellIds = (k-1)*ppf*nWell+1:(nWell*ppf)
    if( k==ceiling(nPlate/ppf))
      wellIds=(k-1)*ppf*nWell+1:(nWell*(nPlate-ppf*(k-1)))
    y.vec = inten.vec0[wellIds]
    if( is.null(yRange) ) yRange = range(y.vec, na.rm=T)  
    theX.df = x.df[wellIds, ]
    plot(range(theX.df$x), yRange, type="n", axes=F,
         xlab=ifelse(is.na(xlab), "", xlab),
         ylab=ifelse(is.na(ylab), intensityName, ylab), 
         main=ifelse(is.na(main),
              paste(intensityName, "-", ifelse(byRow, "By Row", "By Column")),
              main)
         )
    axis(1, at=unique(theX.df[,"plateOrder"])*nWell,
         labels=paste(unique(theX.df[,"plateOrder"]),
                     unique(theX.df[,"plateName"]), sep=": "),
         las=2, cex.axis=0.5 )
    axis(2, las=1); box()
    for( i in order.wellTypes ) {
      condt = theX.df[,"wellType"] == show.wellTypes[i]
      points(theX.df[condt,"x"], y.vec[condt], col=color.wells[i],
             pch=pch.wells[i], cex=cex.point)
  }
  legend( ifelse(is.na(x.legend), -nPlate+(k-1)*ppf*nWell, x.legend),
          ifelse(is.na(y.legend), yRange[2], y.legend), legend=show.wellTypes,
          col=color.wells, pch=pch.wells, cex=cex.legend)
  }
}
