require("data.table")
require("fftwtools")
require("pracma")
require("gstat")
require("sp")
require("stringr")
#require(reshape2)

if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("x", "y"))


#' AFM image class
#' 
#' A S4 class to store and manipulate images from Atomic Force Microscopes.
#'
#' @slot data ($x,$y,$h): a data.table storing the coordinates of the sample and the measured heights
#' @slot samplesperline number of samples per line (e.g.: 512)
#' @slot lines number of line (e.g.: 512)
#' @slot hscansize horizontal size of scan usualy in nanometer (e.g.: hscansize=1000 for a scan size of 1000 nm)
#' @slot vscansize vertical size of scan usualy in nanometer (e.g.: vscansize=1000 for a scan size of 1000 nm)
#' @slot scansize if hscansize equals vscansize, scansize is the size of scan usualy in nanometer (e.g.: scansize=1000 for a scan size of 1000 nm)
#' @slot fullfilename directory and filename on the disk (e.g.: /users/ubuntu/flatten-image.txt)
#' @author M.Beauvais
#' @examples
#' library(AFM)
#' library(data.table)
#' 
#' # create a 128 pixels by 128 pixels AFM image
#' Lines=128
#' Samplesperline=128
#' fullfilename="RandomFakeAFMImage"
#' # the size of scan is 128 nm
#' ScanSize=128
#' # the heights is a normal distribution in nanometers
#' nm<-c(rnorm(128*128, mean=0, sd=1 )) 
#' 
#' scanby<-ScanSize/Samplesperline
#' endScan<-ScanSize*(1-1/Samplesperline)
#' RandomFakeAFMImage<-AFMImage(
#'      data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
#'                        y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
#'                        h = nm),
#'      samplesperline = Samplesperline, lines = Lines, 
#'      vscansize = ScanSize, hscansize = ScanSize, scansize = ScanSize, 
#'      fullfilename = fullfilename )
#' @name AFMImage-class
#' @rdname AFMImage-class
#' @exportClass AFMImage
AFMImage<-setClass("AFMImage",
                   slots = c(data="data.table", 
                             samplesperline="numeric", 
                             lines="numeric", 
                             hscansize="numeric", 
                             vscansize="numeric", 
                             scansize="numeric", 
                             fullfilename="character"
                   ))

#' Constructor method of AFMImage Class.
#'
#' @param .Object an AFMImage object
#' @param data ($x,$y,$h): a data.table storing the coordinates of the sample and the measured heights
#' @param samplesperline number of samples per line (e.g.: 512)
#' @param lines number of line (e.g.: 512)
#' @param hscansize horizontal size of scan usualy in nanometer (e.g.: hscansize=1000 for a scan size of 1000 nm)
#' @param vscansize vertical size of scan usualy in nanometer (e.g.: vscansize=1000 for a scan size of 1000 nm)
#' @param scansize if hscansize equals vscansize, scansize is the size of scan usualy in nanometer (e.g.: scansize=1000 for a scan size of 1000 nm)
#' @param fullfilename directory and filename on the disk (e.g.: /users/ubuntu/flatten-image.txt)
#' @rdname AFMImage-class
#' @export
setMethod(f="initialize",
          signature="AFMImage", 
          definition= function(.Object,
                               data,
                               samplesperline, 
                               lines, 
                               hscansize, 
                               vscansize, 
                               scansize, 
                               fullfilename)  
          { 
            if (!missing(data)) .Object@data<-data
            if (!missing(samplesperline)) .Object@samplesperline<-samplesperline
            if (!missing(lines)) .Object@lines<-lines
            if (!missing(hscansize)) .Object@hscansize<-hscansize
            if (!missing(vscansize)) .Object@vscansize <-vscansize
            if (!missing(scansize)) .Object@scansize <-scansize
            if (!missing(fullfilename)) .Object@fullfilename<-fullfilename
            validObject(.Object)      
            return(.Object)
          })

#' Wrapper function AFMImage
#'
#' @rdname AFMImage-class
#' @export
AFMImage <- function(data,
                     samplesperline, 
                     lines, 
                     hscansize, 
                     vscansize, 
                     scansize, 
                     fullfilename) {
  return(new("AFMImage", data,
             samplesperline, 
             lines, 
             hscansize, 
             vscansize, 
             scansize, 
             fullfilename))
}

#' AFM image sample
#'
#' A real dataset containing an \code{\link{AFMImage}} of an Aluminium interface.
#' The image is made of 512*512 samples of a 1000 nm * 1000 nm surface.
#' samplesperline=512
#' lines=512
#' hscansize=1000
#' vscansize=1000
#' 
#' @name AFMImageOfAluminiumInterface
#' @author J.Landoulsi, I.Liascukiene
NULL

#' AFM image sample
#'
#' A fake dataset containing a manually generated \code{\link{AFMImage}} (peaks regularly positioned on the surface).
#' The image is made of 128*128 samples of a 128 nm * 128 nm surface.
#' samplesperline= 128
#' lines= 128
#' hscansize= 128
#' vscansize= 128
#'
#' @name AFMImageOfRegularPeaks
NULL

#' AFM image sample
#'
#' A fake dataset containing a manually generated \code{\link{AFMImage}} (one peak positioned on the surface).
#' The image is made of 128*128 samples of a 128 nm * 128 nm surface.
#' samplesperline= 128
#' lines= 128
#' hscansize= 128
#' vscansize= 128
#'
#' @name AFMImageOfOnePeak
NULL

#' AFM image sample
#'
#' A fake dataset containing a manually generated \code{\link{AFMImage}} (a normal distribution of heights).
#' The image is made of 128*128 samples of a 128 nm * 128 nm surface.
#' samplesperline= 128
#' lines= 128
#' hscansize= 128
#' vscansize= 128
#'
#' @name AFMImageOfNormallyDistributedHeights
NULL

#' Import data from nanoscope analysis(tm) tool
#'
#' The imported file should contain a header and list of heights
#' The header should contain the following fields:
#' \itemize{
#'   \item Lines: number of scanned lines (e.g. 512)
#'   \item Sampsline: number of scan per line (e.g. 512)
#'   \item ScanSize: the sample size (e.g. 1000nm) the extension nm is mandatory and will be removed
#' }
#' 
#' \code{importFromNanoscope} returns an \code{\link{AFMImage}}
#' @param fullfilename a concatenated string of the directory and filename exported with Nanoscope analysis(TM) software
#' @author M.Beauvais
#' @rdname AFMImage-importFromNanoscope
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' fullfilename<-"/user/ubuntu/NanoscopeFlattenExportedFile.txt"
#' myAFMimage<-importFromNanoscope(fullfilename)
#' displayIn3D(myAFMimage, 1024)
#' }
importFromNanoscope<-function(fullfilename){
  
  #print(fullfilename)
  filename =basename(fullfilename)
  #print(filename)
  # "\Version: 0x08100000"
  # "\Samps/line: 512"
  # "\Lines: 512"
  # "\Scan Size: 1000 nm"
  #Height(nm)
  headerEndString<-"Height(nm)"
  
  wholeFile <- fread(fullfilename)
  wholeFile <- unlist(wholeFile)
  headerSizeWhich<-which(wholeFile==headerEndString)+1
  #print(headerSizeWhich)
  
  hdrs <- read.table(fullfilename, 
                     nrows=headerSizeWhich, 
                     skip=2,
                     comment.char="", 
                     strip.white= TRUE, 
                     check.names=TRUE)
  newhdrs<-sapply(hdrs, function(x) {
    x<-str_replace_all(x, pattern="[^0-9a-zA-Z,.:]+", replacement="")
  })
  
  unlistedNewHdrs=unlist(strsplit(newhdrs, ":", fixed=TRUE))
  
  oneSamplesperline =which(tolower(unlistedNewHdrs)== "sampsline")
  Samplesperline = as.numeric(unlistedNewHdrs[oneSamplesperline+1][1])
  #Samplesperline
  oneScanSize =which(tolower(unlistedNewHdrs)== "scansize")
  ScanSize= unlistedNewHdrs[oneScanSize+1][1]
  ScanSize = as.numeric(substr(ScanSize, 1, nchar(ScanSize)-2))
  #ScanSize
  oneLines =which(tolower(unlistedNewHdrs)== "lines")
  Lines = as.numeric(unlistedNewHdrs[oneLines+1][1])
  #Lines
  
  
  print(paste(ScanSize, Samplesperline, Lines))
  nM <- read.table(fullfilename, skip=headerSizeWhich)
  nM <- unlist(nM)
  
  scanSizeFromZero<-ScanSize-1
  scanby<-ScanSize/Samplesperline
  endScan<-ScanSize*(1-1/Samplesperline)
  
  print(paste("imported ", filename, "...", sep=""))
  
  return(new ("AFMImage", data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
                                            y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
                                            h = nM), 
              samplesperline = Samplesperline, 
              lines = Lines, 
              hscansize = ScanSize,
              vscansize = ScanSize,
              scansize = ScanSize,
              fullfilename = fullfilename))
}


#' Save an AFM image on disk.
#' 
#' The function saves the an \code{\link{AFMImage}} as a rdata file. It uses the fullfilename param of the \code{\link{AFMImage}} and add "AFMImage.rda" extension to save the rdata file on disk.
#' 
#' \code{saveOnDisk} save on disk an \code{\link{AFMImage}} as rdata file
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param exportDirectory an optional argument to change the directory where the rdata file will be stored on disk
#' @author M.Beauvais
#' @rdname AFMImage-saveOnDisk
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' # save the rdata file of the AFMImage in the tempdir() directory;
#' # select another directory to save it permanently on your hard drive
#' saveOnDisk(AFMImageOfAluminiumInterface, tempdir())
#' 
saveOnDisk<-function(AFMImage, exportDirectory){
  if (missing(exportDirectory)) {
    exportDirectory=dirname(AFMImage@fullfilename)
  }
  fullfilename<-paste(exportDirectory, paste(basename(AFMImage@fullfilename),"AFMImage.rda",sep="-"), sep="/")
  save(AFMImage, file=fullfilename)
}

#' Get a sample of an AFM image.
#' 
#' Random selection of heights to keep in an \code{\link{AFMImage}}.
#' This function can be used to calculate quickly an approximated variogram of a large image.
#' 
#' \code{sampleAFMImage} returns a sample of the AFMImage to boost calculation time of variogram
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param percentage percentage of heights to keep
#' @return a sample of an \code{\link{AFMImage}}
#' @author M.Beauvais
#' @rdname AFMImage-sampleAFMImage
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' library(ggplot2)
#' 
#' data(AFMImageOfAluminiumInterface)
#' anAFMImageSample<-sampleAFMImage(AFMImageOfAluminiumInterface,15)
#' variogramAnalysis<-AFMImageVariogramAnalysis(sampleFitPercentage=3.43)
#' avario<-AFM::calculateOmnidirectionalVariogram(AFMImage= anAFMImageSample, 
#'                                                AFMImageVariogramAnalysis= variogramAnalysis)
#' dist<-gamma<-NULL
#' p1 <- ggplot(avario, aes(x=dist, y=gamma))
#' p1 <- p1 + geom_point()
#' p1 <- p1 + geom_line()
#' p1 <- p1 + ylab("semivariance")
#' p1 <- p1 + xlab("distance (nm)")
#' p1 <- p1 + ggtitle("Approximation of variogram thanks to sampling")
#' p1
#' }
#' 
sampleAFMImage<-function(AFMImage, percentage) {
  totalSize<-nrow(AFMImage@data)
  sampleSize<-floor(totalSize*percentage/100)
  
  sampleAFMImage<-AFMImage
  sampleAFMImage@data<-AFMImage@data[sample(1:totalSize, size=sampleSize, replace=F), ]
  sampleAFMImage@fullfilename<-paste(sampleAFMImage@fullfilename, "-sample.txt", sep="")
  sampleAFMImage
}

#' Extract a portion of an AFM image.
#' 
#' The extract will be a square of the specified size.
#' If the size is too large for the original \code{\link{AFMImage}}, only the biggest valid size will be kept.
#' 
#' \code{extractAFMImage} returns an extract of the AFMImage
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param cornerX horizontal coordinates of the extract
#' @param cornerY vertical coordinates of the extract
#' @param size square size of the extract
#' @return a new  \code{\link{AFMImage}} sample
#' @author M.Beauvais
#' @export
#' @rdname AFMImage-extractAFMImage
#' @examples
#' data(AFMImageOfAluminiumInterface)
#' anAFMImageExtract<-extractAFMImage(AFMImageOfAluminiumInterface,15,15,256)
#' 
extractAFMImage<-function(AFMImage, cornerX, cornerY, size) {
  size<-size-1
  
  cornerX<-cornerX*AFMImage@hscansize/AFMImage@samplesperline
  cornerY<-cornerY*AFMImage@vscansize/AFMImage@lines
  size=size*AFMImage@hscansize/AFMImage@samplesperline
  
  minX<-cornerX
  maxX<-cornerX+size
  
  minY<-cornerY
  maxY<-cornerY+size
  print(paste(cornerX,cornerY,size,minX, maxX, minY,maxY))
  
  alldata<-copy(AFMImage@data)
  key(alldata)
  keycols = c("y","x")
  setkeyv(alldata,keycols) 
  
  x<-y<-NULL
  alldata<-alldata[x>=minX & x<=maxX & y>=minY & y<=maxY]
  alldata$x<-alldata$x-min(alldata$x)
  alldata$y<-alldata$y-min(alldata$y)
  vscansize<-max(alldata$x)-min(alldata$x)+1
  hscansize<-max(alldata$y)-min(alldata$y)+1
  scansize<-max(vscansize, hscansize)
  samplesperline<-length(unique(alldata$y))
  lines<-length(unique(alldata$x))
  fullfilename<-paste(AFMImage@fullfilename, "extract.txt", sep="-")
  
  
  newAFMImage=new("AFMImage",
                  data=alldata,
                  vscansize=vscansize,
                  hscansize=hscansize,
                  scansize=scansize,
                  samplesperline=samplesperline,
                  lines=lines,
                  fullfilename=fullfilename
  )
  
  
  #   newAFMImage@data<-copy(AFMImage@data)
  #   alldata<-newAFMImage@data
  #   key(alldata)
  #   keycols = c("y","x")
  #   setkeyv(alldata,keycols) 
  #   #max(AFMImage@data$x)
  #   #max(alldata$h)
  #   newAFMImage@data<-alldata[x>=minX & x<=maxX & y>=minY & y<=maxY]
  #   newAFMImage@data$x<-newAFMImage@data$x-min(newAFMImage@data$x)
  #   newAFMImage@data$y<-newAFMImage@data$y-min(newAFMImage@data$y)
  #   newAFMImage@vscansize<-max(newAFMImage@data$x)-min(newAFMImage@data$x)+1
  #   newAFMImage@hscansize<-max(newAFMImage@data$y)-min(newAFMImage@data$y)+1
  #   newAFMImage@scansize<-max(newAFMImage@vscansize, newAFMImage@hscansize)
  #   newAFMImage@samplesperline<-length(unique(newAFMImage@data$y))
  #   newAFMImage@lines<-length(unique(newAFMImage@data$x))
  #   newAFMImage@fullfilename<-paste(AFMImage@fullfilename, "extract.txt", sep="-")
  return(newAFMImage)
}


#' simplify an AFM image.
#' 
#' The simplification is taking a very simple gridded sample of the image. 
#' It can be useful to speed up display.
#' 
#' \code{simplifyAFMImage} returns a simplified AFMImage
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param newSamplesperline the new number of samplesperline of the AFMImage
#' @param newLines the new number of lines of the AFMImage
#' @return a new simplified \code{\link{AFMImage}}
#' @author M.Beauvais
#' @export
#' @rdname AFMImage-simplifyAFMImage
#' @examples
#' data(AFMImageOfAluminiumInterface)
#' anAFMImageExtract<-simplifyAFMImage(AFMImageOfAluminiumInterface,16,16)
#' 
simplifyAFMImage<-function(AFMImage, newSamplesperline, newLines) {
  print(paste("simplifyAFMImage", newSamplesperline, newLines))

  #   AFMImage<-AFMImageOfRegularPeaks
  #   newSamplesperline<-16
  #   newLines<-16
  # print(paste(newSamplesperline, newLines))
  
  if (newSamplesperline> AFMImage@samplesperline) newSamplesperline<-AFMImage@samplesperline
  if (newLines> AFMImage@lines) newLines<-AFMImage@lines
  
  
  z<-matrix(AFMImage@data$h,nrow = AFMImage@lines,ncol = AFMImage@samplesperline)
  samplesperlineBy=ceil(AFMImage@samplesperline/newSamplesperline)
  samplesperlineIndices=seq(1,AFMImage@samplesperline, by=samplesperlineBy)
  
  linesBy=ceil(AFMImage@lines/newLines)
  linesIndices=seq(1,AFMImage@lines, by=linesBy)
  newZ=z[samplesperlineIndices,linesIndices]
  
  
#   print(paste(samplesperlineBy))
#   print(paste(samplesperlineIndices))
#   print(paste(linesBy))
#   print(paste(linesIndices))
#   print(newZ)
  
  
  scanby<-AFMImage@scansize/newSamplesperline
  endScan<-AFMImage@scansize*(1-1/newSamplesperline)
  
  newData = data.table(x = rep(seq(0,endScan, by= scanby), times = newLines),
                       y = rep(seq(0,endScan, by= scanby), each = newSamplesperline),
                       h =as.numeric(newZ))
  
#   print(newData$x)
#   print(newData$y)
  
  newAFMImage=new("AFMImage",
                  data=newData,
                  vscansize=AFMImage@vscansize,
                  hscansize=AFMImage@hscansize,
                  scansize=AFMImage@scansize,
                  samplesperline=newSamplesperline,
                  lines=newLines,
                  fullfilename=AFMImage@fullfilename)
  
  return(newAFMImage)
  
}
