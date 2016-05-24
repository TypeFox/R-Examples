imageIntensity.fn <-
function(data.df, intensityName=NA, plateName="BARCODE", wellName="WELL_USAGE",
         rowName, colName, sampleName="Sample", sourcePlateName=NA)
{
#*************************************************************************#
# author: Xiaohua Douglas Zhang & Zhaozhi Zhang, 2012                     #
# function to display the image of intensity of each well in a plate      #
# Input:                                                                  #
#   data.df: data for all plates including at least five columns          #
#       for well names, rows and columns, respectively                    #
#   wellName: data column name indicating well names                      #
#   rowName:  data column name indicating row numbers in a plate          #
#   colName:  data column name indicating column numbers in a plate       #
#   wells:  names for unique wells                                        #
#   colors:  colors for corresponding unique wells                        #
#   sourcePlateName: column name in data.df indicating soure plate names  #
# References:
#   Zhang XHD, Zhang ZZ. 2012. displayHTS: a R package displaying data and
#        results from high-throughput screening experiments (submitted).  
#   Zhang XHD, 2011. Optimal High-Throughput Screening: Practical Experimental
#        Design and Data Analysis for Genome-scale RNAi Research.
#        Cambridge University Press, Cambridge, UK.
#   Zhang XHD, 2008. Novel analytic criteria and effective plate designs for
#        quality control in genome-wide RNAi screens.
#        Journal of Biomolecular Screening 13(5): 363-377.
#   Zhang XHD, Espeseth AS, Johnson E, Chin J, Gates A, Mitnaul L, Marine SD,
#        Tian J, Stec EM, Kunapuli P, Holder DJ, Heyse JF, Stulovici B,
#        Ferrer M. 2008. Integrating experimental and analytic approaches to
#        improve data quality in genome-wide RNAi screens.
#        Journal of Biomolecular Screening 13(5): 378-389.   
# See Also:
#   imageDesign.fn, dualFlashlight.fn, plateWellSeries.fn
# Example:                                                                #
#   data("HTSdataSort", package = "displayHTS")
#   par( mfrow=c(4, 3) )
#   imageIntensity.fn(HTSdataSort, intensityName="log2Intensity",             
#                      plateName="BARCODE", wellName="WELL_USAGE",
#                      rowName="XPOS", colName="YPOS", 
#                      sampleName="Sample", sourcePlateName="SOBARCODE" ) 
#*************************************************************************#

  fac2char.fn <- function(x) {levels(x)[as.numeric(x)]}

  getIntensityMatrix.fn <- function( dataIn.df, intensity=T, nrow=16, ncol=24)
  {
  ##----------------------------------------------------------------
  ## Functions to get the intensity matrix
  ## dataIn.df must have three columns respectively corresponding to
  ##  row position index, column position index and response
  ## The plate must have 'nrow' rows and 'ncol' columns
  # author: Xiaohua Douglas Zhang, 2005                                  
  ##-----------------------------------------------------------------
      if( dim(dataIn.df)[1] != nrow*ncol ) { 
        new.mat = matrix(NA, nrow, ncol)
        if( intensity ) {
          response.vec = dataIn.df[,3]
        } else {
          response.vec = fac2char.fn(dataIn.df[,3])
        }
        for(i in 1:nrow)
          for(j in 1:ncol) {
            idx = dataIn.df[,1]==i& dataIn.df[,2]==j 
            if( sum(idx) == 1 ) 
              new.mat[i,j] = response.vec[idx]
          }
      } else {
        data.df = dataIn.df[ order(dataIn.df[,1], dataIn.df[,2]),]
        new.mat = matrix(data.df[,3], nrow, ncol, byrow=T)
      }
      new.mat
  }

  matrixImage.fn <-
  function(matr, zlim.low=min(matr, na.rm=T), zlim.high=max(matr, na.rm=T),
           col=rgb(red=(0:20)/20, green=(20:0)/20, blue=0), xlab="", ylab="",
           cex.axis=1, gridcol="grey", gridlwd=1, cex.outlier=2,
           colUp="red", colDown="green", labUp="+", labDown="-")
  {
    x=1:ncol(matr)
    y=1:nrow(matr)
    zlim=c(zlim.low, zlim.high)
    
    inten.mat = t(matr)[,nrow(matr):1]
    upI = inten.mat > zlim.high
    downI = inten.mat < zlim.low
    if(sum( upI, na.rm=T) != 0 ) inten.mat[upI] = NA
    if(sum( downI, na.rm=T) != 0 ) inten.mat[downI] = NA
    image(x, y, inten.mat, zlim=zlim, col=col, axes=FALSE, xlab=xlab, ylab=ylab)
    if(sum( upI, na.rm=T) != 0 )
      text( row(inten.mat)[upI], col(inten.mat)[upI], labels=labUp,
           col=colUp, cex=cex.outlier)
    if(sum( downI, na.rm=T) != 0 )  
      text( row(inten.mat)[downI], col(inten.mat)[downI], labels=labDown,
           col=colDown, cex=cex.outlier)
    grid(ncol(matr),nrow(matr),col=gridcol,lty = 1,lwd=gridlwd) 
    axis(2,at=1:nrow(matr),labels=nrow(matr):1,lwd=0.5,las=1,cex.axis=cex.axis)
    axis(3,at=1:ncol(matr),labels=1:ncol(matr),lwd=0.5,las=2,cex.axis=cex.axis)
    box()
  }

  image.scale <-  
  function(sample=1, z, col, x, y = NULL, size = NULL, digits = 2,
           labels = c("breaks", "ranges"), bkg.names, height.factor=0.618,
           cex.label=1)
   {
     n = length(col)
     usr = par("usr")
     mx = mean(usr[1:2]); my = mean(usr[3:4])
     dx = diff(usr[1:2]); dy = diff(usr[3:4])
     if (missing(x))
       x = mx + 1.05*dx/2 
     else if (is.list(x)) {
       if (length(x$x) == 2) 
         size = c(diff(x$x), -diff(x$y)/n)
       y = x$y[1]
       x = x$x[1]
     } else x = x[1]
     if (is.null(size))
       if (is.null(y)) {
         size = height.factor*dy/n   
         y = my + height.factor*dy/2  
       } else size = (y-my)*2/n
     if (length(size)==1)
       size = rep(size, 2)    
     if (is.null(y))
       y = my + n*size[2]/2
     
     i = seq(along = col)
     rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
          col = rev(col), xpd = TRUE)
     
     rng = range(z, na.rm = TRUE)
     bks = seq(from = rng[2], to = rng[1], length = n + 1)
     bks = formatC(bks, format="f", digits=digits)
     labels = match.arg(labels)
     if (labels == "breaks")
       ypts = y - c(0, i) * size[2]
     else {
       bks = paste(bks[-1], bks[-(n+1)], sep = " - ")
       ypts = y - (i - 0.5) * size[2]
     }
    
     if (sample == 1){
       text(x = x + 1.2 * size[1], y = ypts, labels = bks,
            adj = ifelse(size[1]>0, 0, 1), xpd = TRUE, cex=cex.label) 
     }
     else{
       bks = bkg.names
       text(x = x + 1.2 * size[1], y = ypts+size[2]/2, labels = bks,
            adj = ifelse(size[1]>0, 0, 1), xpd = TRUE, cex=cex.label) 
     }
    
  }

  imageContinuous.fn <- 
  function(matr, continuous=1,
           zlim.low=min(matr, na.rm=T),zlim.high=max(matr, na.rm=T),
           col=rgb(red=(0:20)/20, green=(20:0)/20, blue=0), 
           gridcol="lightgray", gridlwd=1,
           cex.outlier=2, colUp="red", colDown="green", labUp="+", labDown="-",
           mar=c(4, 4, 5, 5) + 0.1, height.factor=0.9, cex.label=1, 
           xlab="", ylab="", cex.axis=1 )
  {
    par("mar" = mar)       
    matrixImage.fn(matr,zlim.low=zlim.low, zlim.high=zlim.high,
                   col=col,xlab=xlab, ylab=ylab, cex.axis=cex.axis,
                   gridcol=gridcol, gridlwd=gridlwd, cex.outlier=cex.outlier,
                   colUp=colUp, colDown=colDown, labUp=labUp, labDown=labDown)
    matr1 = matr
    matr1[matr <= zlim.low ] = zlim.low
    matr1[matr >= zlim.high] = zlim.high  
    image.scale(1,matr1,col=col,height.factor=height.factor,cex.label=cex.label)
    par(mar=c(5.1, 4.1, 4.1, 2.1)) 
  }

  inten.vec0= data.df[, intensityName]
  plate.vec = as.vector(data.df[, plateName])
  plates = unique( plate.vec )
  nPlate = length(plates)
  if(!is.na(sourcePlateName)) {
    soPlate.vec = as.vector(data.df[, sourcePlateName])
    platesBar=unique(cbind(soPlate.vec, plate.vec))
  }

  nRow=max(data.df[,rowName]); nCol=max(data.df[,colName])

  for( i in 1:nPlate ) { 
    condt = plate.vec==plates[i]
    dataIn.df = data.df[ condt, c(rowName,colName,wellName) ]
    wellusage = dataIn.df[, wellName]
    inten.vec = inten.vec0[ condt ]
    sampleCondt = wellusage=="Sample"
    boxplot.boundary =
      boxplot(inten.vec[sampleCondt], plot=F)$stats[c(1,5)]
    inten.mat =
      getIntensityMatrix.fn(cbind(dataIn.df[, c(rowName, colName)], inten.vec),
                            nrow=nRow, ncol=nCol)
    imageContinuous.fn(inten.mat, zlim.low = boxplot.boundary[1],
                        zlim.high = boxplot.boundary[2], gridlwd=0.25,
                        mar=c(1, 1.5, 2, 1.5) + 0.1, height.factor=1.1,
                        cex.outlier=0.5, cex.label=0.3, cex.axis=0.5)
    if( is.na(sourcePlateName) ) {
      title( sub=plates[i] )
    } else {
      title( sub=paste(platesBar[i,1], "-", platesBar[i,2], "         ") )
    }
  }
}
