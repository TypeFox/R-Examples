require("fftwtools")
require("pracma")

require("data.table")

require("gstat")
require(sp)

require(rgl)
require(reshape2)

setOldClass("mesh3d")

#' @title AFM image Power Spectrum Density analysis class
#' 
#' @description \code{AFMImage3DModelAnalysis} 
#'
#' @slot f1 a face of the 3D model
#' @slot f2 a face of the 3D model
#' @slot f3 a face of the 3D model
#' @slot f4 a face of the 3D model
#' @name AFMImage3DModelAnalysis-class
#' @rdname AFMImage3DModelAnalysis-class
#' @author M.Beauvais
AFMImage3DModelAnalysis<-setClass("AFMImage3DModelAnalysis",
                                  slots = c(
                                    f1="mesh3d",
                                    f2="mesh3d",
                                    f3="mesh3d",
                                    f4="mesh3d",
                                    updateProgress="function"),
                                  validity = function(object) { 
                                    return(TRUE)
                                  })

#' Display a 3D image of an AFMImage and store it on disk.
#'
#' Display a 3D image of an AFMImage and store it on disk if fullfilename variable is set.
#' It uses the \code{\link{rgl}} package.
#' 
#' @param AFMImage the AFM image to be displayed in three dimensions.
#' @param fullfilename (optional) the directory and filename to save the png of the 3D image. If this variable is missing, the function will not save on disk the 3D image.
#' @param width (optional)  width of the image. Default is 512 pixels. Note: width can't be superior to screen resolution.
#' @param changeViewpoint (optional) if TRUE, the viewpoint is changed. Default is TRUE.
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' newAFMImage<-extractAFMImage(AFMImageOfAluminiumInterface, cornerX=50, cornerY=50, size=256)
#' displayIn3D(newAFMImage, 1024)
#'
#' \dontrun{
#' AFMImage<-importFromNanoscope("/user/ubuntu/myNanoscopeAnalysisAFMImage.txt")
#' displayIn3D(AFMImage, 1024, "/users/ubuntu/myRglAFMimage.png")
#' }
#' 
displayIn3D<- function(AFMImage, width, fullfilename, changeViewpoint) {
  if(missing(width)){
    width <- 512
  }
  if(missing(fullfilename)){
    save <- FALSE
  }else{
    save <- TRUE
  }
  if(missing(changeViewpoint)){
    changeViewpoint<-TRUE
  }else{
    changeViewpoint<-FALSE
  }
  
  # respect the proportion between horizontal / vertical distance and heigth
  newHeights <- (AFMImage@data$h)*(AFMImage@samplesperline)/(AFMImage@scansize)
  
  minH<-min(newHeights)
  # TODO check validity of created image instead
  if(!is.na(minH)) {
    newH<-(newHeights-minH)
    y<-matrix(newH, nrow = AFMImage@lines, ncol = AFMImage@samplesperline)
    #z <- seq(ncol(y),1,by=-1) 
    z <- seq(1,ncol(y),by=1) 
    x <- (1:nrow(y))
    
    ylim <- range(y)
    ylen <- ylim[2] - ylim[1] + 1
    print(ylen)
    colorlut <- heat.colors(ylen, alpha = 1) # height color lookup table
    col <- colorlut[ y-ylim[1]+1 ] # assign colors to heights
    
    rgl.open()
    par3d(windowRect = 100 + c( 0, 0, width, width ) )
    rgl.bg(color = c("white"),  back = "lines")
    
    bboxylen=3
    if(ylim[2]<60) bboxylen=2
    
    rgl.bbox(color = c("#333333", "black"), emission = "#333333", 
             specular = "#111111", shininess = 0, alpha = 0.6, xlen=0, zlen=0, ylen=bboxylen )
    rgl.surface(x, z, y, color=col, back="lines")
    
    if (changeViewpoint) {
      i<-130
      rgl.viewpoint(i,i/4,zoom=1.1)
    }
    if (save) {
      print(paste("saving", basename(fullfilename)))
      rgl.snapshot(fullfilename)
    }
    return(TRUE)
  }
  return(FALSE)
}



#' Calculate the 3D model for 3D printing
#'
#' \code{calculate3DModel} update  \code{\link{AFMImage3DModelAnalysis}}
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImage3DModelAnalysis n \code{\link{AFMImage3DModelAnalysis}} to store the setup and results of PSD analysis
#' 
#' @name calculate3DModel
#' @rdname calculate3DModel-methods
#' @exportMethod calculate3DModel
#' @author M.Beauvais
setGeneric(name= "calculate3DModel", 
           def= function(AFMImage3DModelAnalysis, AFMImage) {
             return(standardGeneric("calculate3DModel"))
           })

#' @rdname calculate3DModel-methods
#' @aliases calculate3DModel,AFMImage-method
setMethod(f="calculate3DModel", "AFMImage3DModelAnalysis",
          definition= function(AFMImage3DModelAnalysis, AFMImage) {

            print(paste("exporting to stl format "))
            baseThickness<-2
            
            # respect the proportion between horizontal / vertical distance and heigth
            newHeights <- (AFMImage@data$h)*(AFMImage@samplesperline)/(AFMImage@scansize)
            
            minH<-min(newHeights)
            #print(paste("minH", minH))
            if (minH<0) { newH<-(newHeights-minH+baseThickness) 
            }  else { newH<-(newHeights-minH+5) }
            #print(paste("min(newH)", min(newH)))
            #print(paste("max(newH)", max(newH)))
            
            totalLength<-4
            counter<-0
            
            
            if (!is.null(AFMImage3DModelAnalysis@updateProgress)&&
                is.function(AFMImage3DModelAnalysis@updateProgress)&&
                !is.null(AFMImage3DModelAnalysis@updateProgress())) {
              text <- paste0("starting ", totalLength, " calculations")
              #AFMImage3DModelAnalysis@updateProgress(message="Calculating 3D faces", value=0)
              
              AFMImage3DModelAnalysis@updateProgress(value= 0, detail = text)
              
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0(round(counter, 2),"/",totalLength)
              AFMImage3DModelAnalysis@updateProgress(value= value, detail = text)
              print("update")
            }else{
              print("no GUI update")
              print(is.null(AFMImage3DModelAnalysis@updateProgress))
              print(is.function(AFMImage3DModelAnalysis@updateProgress))
              print(is.null(AFMImage3DModelAnalysis@updateProgress()))
            }
            
            
            #face 1
            x1<-seq(1:AFMImage@lines)
            y1<-rep(rep(1, each = AFMImage@lines) , each=1)
            z1<-newH[x1+(y1-1)*AFMImage @samplesperline]
            
            x1=c(x1,x1[length(x1)])
            y1=c(y1,1)
            z1=c(z1,1)
            
            x1=c(x1,1)
            y1=c(y1,1)
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,z1[1])
            
            #    print(length(x1))
            #    print(length(y1))
            #    print(length(z1))
            #    print(x1)
            #    print(y1)
            #    print(z1)
            
            f1<-polygon3d(x1, z1, y1, col = "red", plot=FALSE, fill=TRUE)
            f1<-rotate3d( f1 , -pi/2, 1, 0, 0 )
            f1<-translate3d( f1 , 0, AFMImage@samplesperline+1, 0 )
            
            #face 2
            if (!is.null(AFMImage3DModelAnalysis@updateProgress)&&
                is.function(AFMImage3DModelAnalysis@updateProgress)&&
                !is.null(AFMImage3DModelAnalysis@updateProgress())) {
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0(round(counter, 2),"/",totalLength)
              AFMImage3DModelAnalysis@updateProgress(value= value, detail = text)
            }
            
            y1<-rep(AFMImage@lines, each = AFMImage@lines)
            x1<-seq(1:AFMImage@lines)
            z1<-newH[x1+(y1-1)*AFMImage@samplesperline]
            
            x1=c(x1,x1[length(x1)])
            y1=c(y1,y1[1])
            z1=c(z1,1)
            
            x1=c(x1,1)
            y1=c(y1,y1[1])
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,z1[1])
            
            #   print(length(x1))
            #   print(length(y1))
            #   print(length(z1))
            #   print(x1)
            #   print(y1)
            #   print(z1)
            y1<-as.numeric(y1)
            
            # z1<-z1+rnorm(1:length(z1))
            
            f2<-polygon3d(x1, z1, y1, col = "blue", plot=FALSE, fill=TRUE)
            f2<-rotate3d( f2 , -pi/2, 1, 0, 0 )
            f2<-translate3d( f2 , 0, AFMImage@lines+1, 0 )
            
            
            #face 3
            if (!is.null(AFMImage3DModelAnalysis@updateProgress)&&
                is.function(AFMImage3DModelAnalysis@updateProgress)&&
                !is.null(AFMImage3DModelAnalysis@updateProgress())) {
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0(round(counter, 2),"/",totalLength)
              AFMImage3DModelAnalysis@updateProgress(value= value, detail = text)
            }
            
            y1<-seq(1:AFMImage@samplesperline)
            x1<-rep(1, times = (AFMImage@samplesperline))
            z1<-rev(newH[x1+(y1-1)*AFMImage@samplesperline])
            
            x1=c(x1,x1[length(x1)])
            y1=c(y1,y1[length(y1)])
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,z1[1])
            
            #   print(x1)
            #   print(y1)
            #   print(z1)
            
            f3<-polygon3d(y1, z1, x1, col = "red", plot=FALSE, fill=TRUE)
            f3<-rotate3d( f3 , -pi/2, 1, 0, 0 )
            f3<-rotate3d( f3 , -pi/2, 0, 0, 1 )
            f3<-translate3d( f3 , 0, 0, 0 )
            
            #face 4
            if (!is.null(AFMImage3DModelAnalysis@updateProgress)&&
                is.function(AFMImage3DModelAnalysis@updateProgress)&&
                !is.null(AFMImage3DModelAnalysis@updateProgress())) {
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0(round(counter, 2),"/",totalLength)
              AFMImage3DModelAnalysis@updateProgress(value= value, detail = text)
            }
            
            y1<-seq(1:AFMImage@samplesperline)
            x1<-rep(AFMImage@lines, times = AFMImage@samplesperline)
            z1<-rev(newH[x1+(y1-1)*AFMImage@samplesperline])
            
            
            x1=c(x1,x1[length(x1)])
            y1=c(y1,y1[length(y1)])
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,1)
            
            x1=c(x1,x1[1])
            y1=c(y1,y1[1])
            z1=c(z1,z1[1])
            
            f4<-polygon3d(y1, z1, x1, col = "red", plot=FALSE, fill=TRUE)
            f4<-rotate3d( f4 , -pi/2, 1, 0, 0 )
            f4<-rotate3d( f4 , -pi/2, 0, 0, 1 )
            
            AFMImage3DModelAnalysis@f1<-f1
            AFMImage3DModelAnalysis@f2<-f2
            AFMImage3DModelAnalysis@f3<-f3
            AFMImage3DModelAnalysis@f4<-f4
            return(AFMImage3DModelAnalysis)
          })



#' Export an AFM Image as a STL format file.
#' 
#' Export an \code{\link{AFMImage}} as a STL format file thanks to the \code{\link{rgl}} package. The STL file can be used as an input for a 3D printing software tool.\cr\cr
#' exportToSTL is compatible with slicr (http://slic3r.org) version 1.2.9 (GPL v3 licence).\cr
#' In order to 3D print the AFM Image with slic3r, do as following:
#' \itemize{
#' \item Use "File> Repair STL file..." menu option to create a file with the obj extension.
#' \item Use "Add" button below the menu to display your AFM Image on the print board
#' \item Right click on your AFM image. Use "Scale> uniformely" option, Set "15%"  for your AFM image to fit your printing board
#' }
#' @param AFMImage3DModelAnalysis an \code{\link{AFMImage3DModelAnalysis}}
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param stlfullfilename  directory and filename to save as a stl file
#' @author M.Beauvais
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' newAFMImage<-extractAFMImage(AFMImageOfAluminiumInterface, cornerX=50, cornerY=50, size=64)
#' exportToSTL(newAFMImage, paste(tempdir(), "myFile.stl", sep="/"))
#' }
exportToSTL<- function(AFMImage3DModelAnalysis, AFMImage, stlfullfilename) {
  
  print(paste("exporting to stl format ", basename(stlfullfilename) ))
  baseThickness<-2
  
  
  
  #AFMImage3DModelAnalysis<-calculate3DModel(AFMImage3DModelAnalysis= AFMImage3DModelAnalysis, AFMImage= AFMImage)
  
  # respect the proportion between horizontal / vertical distance and heigth
  newHeights <- (AFMImage@data$h)*(AFMImage@samplesperline)/(AFMImage@scansize)
  
  minH<-min(newHeights)
  #print(paste("minH", minH))
  if (minH<0) { newH<-(newHeights-minH+baseThickness) 
  }  else { newH<-(newHeights-minH+5) }
  #print(paste("min(newH)", min(newH)))
  #print(paste("max(newH)", max(newH)))
  
  
  
  # surface
  z<-matrix(newH,nrow = AFMImage@lines,ncol = AFMImage@samplesperline)
  x <- (1:nrow(z)) 
  y <- seq(ncol(z),1,by=-1)
  
  zlim <- range(z)
  zlen <- zlim[2] - zlim[1] + 1
  colorlut <- heat.colors(zlen) 
  col <- colorlut[ z-zlim[1]+1 ]
  
  rgl.open()
  par3d(windowRect = c(100,100,800,800))
  shade3d(AFMImage3DModelAnalysis@f1)
  shade3d(AFMImage3DModelAnalysis@f2)
  shade3d(AFMImage3DModelAnalysis@f3)
  shade3d(AFMImage3DModelAnalysis@f4)
  terrain3d(x, y, z, color=col, front="lines", back="lines")
  
  # create a stl file
  print(paste("saving", basename(stlfullfilename)))
  writeSTL(stlfullfilename)
  print("done")
}

getHorizontalSlice<-function(AFMImage, levelMin, levelMax, width, fullfilename) {
  if(missing(width)){
    width <- 512
  }
  if(missing(fullfilename)){
    save <- FALSE
  }else{
    save <- TRUE
    print(fullfilename)
  }
  
  print("getHorizontalSlice")
  print(paste("levelMin=", levelMin, "levelMax= ",levelMax))
  
  heights<-AFMImage@data$h
  minH<-min(heights)
  #  print(tail(heights, n=10))
  print(paste("min=", min(heights), "max= ",max(heights), "mean=", mean(heights)))
  
  
  indexfmin<-which( heights < levelMin | heights > levelMax)
  print(head(indexfmin, n=10))
  print(length(indexfmin))
  
  
  
  
  minH<-min(AFMImage@data$h)
  newH<-(AFMImage@data$h-minH)
  y<-matrix(newH,nrow = AFMImage@lines,ncol = AFMImage@samplesperline)
  x <- (1:nrow(y)) 
  z <- (1:ncol(y))
  
  ylim <- range(y)
  ylen <- ylim[2] - ylim[1] + 1
  colorlut <- heat.colors(ylen) # height color lookup table
  col <- colorlut[ y-ylim[1]+1 ] # assign colors to heights
  
  oldMinH<-min(newH)
  AFMImage@data$h[indexfmin]=NA
  print("hhhhh")
  print(head(indexfmin, n=1))
  AFMImage@data$h[head(indexfmin, n=1)]=minH
  
  newH<-AFMImage@data$h
  y<-matrix(newH,nrow = AFMImage@lines,ncol = AFMImage@samplesperline)
  
  rgl.open()
  par3d(windowRect = 100 + c( 0, 0, width, width ) )
  rgl.bg(color = c("white"),  back = "lines")
  
  bboxylen=3
  if(ylim[2]<60) bboxylen=2
  
  rgl.bbox(color = c("#333333", "black"), emission = "#333333", 
           specular = "#111111", shininess = 0, alpha = 0.6, xlen=0, zlen=0, ylen=bboxylen )
  rgl.surface(x, z, y, color=col, back="lines")
  i<-130
  rgl.viewpoint(i,i/4,zoom=1.1)
  
  if (save) {
    print(paste("saving", basename(fullfilename)))
    rgl.snapshot(fullfilename)
  }
  
  
  
}


saveHorizontalSlices<-function(AFMImage, numberOfSlices, width, fullfilename) {
  
  heights<-AFMImage@data$h
  print(paste("min=", min(heights), "max= ",max(heights), "mean=", mean(heights)))
  
  minH=min(heights)
  maxH=max(heights)
  sliceHeight = ceil((maxH- minH)/numberOfSlices)
  
  sliceIndex = 0
  for(i in seq(floor(minH), ceil(maxH), by=sliceHeight)) {
    
    sliceIndex<-sliceIndex+1
    sliceName<-sliceIndex
    if (numberOfSlices>9) {
      if (sliceIndex<10) {
        sliceName<-paste("0", sliceIndex, sep="")
      }
    }
    
    newfullfilename = paste(fullfilename, "3D-horizontal-slice",sliceName, "png", sep=".")
    print(newfullfilename)
    
    getHorizontalSlice(AFMImage, i, i + sliceHeight, width, newfullfilename )
    #getHorizontalSlice(AFMImage, i, i + sliceHeight, width)
  }
  
  
}


get3DImageFullfilename<-function(exportDirectory, imageName) {
  fullfilename<-paste(exportDirectory, paste(imageName,"3D.png",sep="."), sep="/")
  return(fullfilename)
}

