require("fftwtools")
require("pracma")

require("data.table")

require("gstat")
require(sp)

require("stringr")

require(fractaldim)

require(reshape2)


#' @title AFM image fractal dimension method class
#' 
#' @description \code{AFMImageFractalDimensionMethod} stores calculation from one fractal dimension method
#'
#' @slot fd_method Two dimensional function names used to evaluate the fractal dimension and fractal scale
#' @slot fd the value of the fractal dimension
#' @slot fd_scale the value of the fractal scale
#' @name AFMImageFractalDimensionMethod-class
#' @rdname AFMImageFractalDimensionMethod-class
#' @author M.Beauvais
#' @seealso \code{\link{fractaldim}}
AFMImageFractalDimensionMethod<-setClass("AFMImageFractalDimensionMethod",
                                         slots = c(fd_method="character", 
                                                   fd="numeric",
                                                   fd_scale="numeric"))

#' Constructor method of AFMImageFractalDimensionMethod Class.
#'
#' @param .Object an AFMImageFractalDimensionMethod object
#' @param fd_method Two dimensional function names used to evaluate the fractal dimension and fractal scale
#' @param fd the value of the fractal dimension
#' @param fd_scale the value of the fractal scale
#' @rdname AFMImageFractalDimensionMethod-class
#' @export
setMethod(f= "initialize", 
          signature= "AFMImageFractalDimensionMethod", 
          definition= function(.Object, fd_method, fd, fd_scale)
          {
            .Object@fd_method<-fd_method
            .Object@fd<-fd
            .Object@fd_scale<-fd_scale
            validObject(.Object)      
            return(.Object)
          })

#' Wrapper function AFMImageFractalDimensionMethod
#' 
#' @rdname AFMImageFractalDimensionMethod-class
#' @export
AFMImageFractalDimensionMethod <- function(fd_method, fd, fd_scale) {
  return(new("AFMImageFractalDimensionMethod", fd_method=fd_method, fd=fd, fd_scale=fd_scale))
}


#' AFM image fractal dimensions analysis class
#' 
#' A S4 class to handle the fractal dimension calculation with several fractal dimension methods
#'
#' @slot fractalDimensionMethods a list of \code{\link{AFMImageFractalDimensionMethod}}
#' @slot csvFullfilename To be removed ?
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImageFractalDimensionsAnalysis-class
#' @rdname AFMImageFractalDimensionsAnalysis-class
#' @exportClass AFMImageFractalDimensionsAnalysis
#' @author M.Beauvais
#'
AFMImageFractalDimensionsAnalysis<-setClass("AFMImageFractalDimensionsAnalysis",
                                            slots = c(fractalDimensionMethods="list", 
                                                      csvFullfilename="character",
                                                      updateProgress="function"))

#' Constructor method of AFMImageFractalDimensionsAnalysis Class.
#'
#' @param .Object an AFMImageFractalDimensionsAnalysis Class
#' @param fractalDimensionMethods a list of \code{\link{AFMImageFractalDimensionMethod}}
#' @param csvFullfilename To be removed ?
#' @rdname AFMImageFractalDimensionsAnalysis-class
#' @export
setMethod("initialize", "AFMImageFractalDimensionsAnalysis", function(.Object, fractalDimensionMethods, csvFullfilename)  
{
  if(!missing(fractalDimensionMethods)) .Object@fractalDimensionMethods<-fractalDimensionMethods
  if(!missing(csvFullfilename)) .Object@csvFullfilename<-csvFullfilename
  validObject(.Object)      
  return(.Object)
})


#' Wrapper function AFMImageFractalDimensionsAnalysis
#'
#' @rdname AFMImageFractalDimensionsAnalysis-class
#' @export
AFMImageFractalDimensionsAnalysis <- function() {
  return(new("AFMImageFractalDimensionsAnalysis"))
}

#' Method \code{fractalDimensionMethods} returns a list of FractalDimensionMethod objects
#' @name AFMImageFractalDimensionsAnalysis-class
#' @rdname AFMImageFractalDimensionsAnalysis-class
#' 
setGeneric("fractalDimensionMethods",function(object){standardGeneric("fractalDimensionMethods")})
setGeneric(name= "fractalDimensionMethods<-", 
           def= function(AFMImageFractalDimensionsAnalysis, value) {
             return(standardGeneric("fractalDimensionMethods<-"))
           })

#' @rdname AFMImageFractalDimensionsAnalysis-class
#' @aliases fractalDimensionMethods
#' @param object a \code{\link{AFMImageFractalDimensionsAnalysis}}
setMethod("fractalDimensionMethods",signature=signature(object='AFMImageFractalDimensionsAnalysis'),
          function(object) {
            return(object@fractalDimensionMethods)
          }
)
setReplaceMethod(f="fractalDimensionMethods", 
                 signature(AFMImageFractalDimensionsAnalysis = "AFMImageFractalDimensionsAnalysis", value = "list"),
                 definition= function(AFMImageFractalDimensionsAnalysis, value) {
                   AFMImageFractalDimensionsAnalysis@fractalDimensionMethods <- value
                   return(AFMImageFractalDimensionsAnalysis)
                 })

#' Calculate 2D fractal dimensions and scales of an AFM Image
#'
#' \code{getFractalDimensions} calculates fractal dimensions and scales of an \code{\link{AFMImage}} with \code{\link[fractaldim]{fd.estim.method}} from the \code{\link{fractaldim}} package.
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageFractalDimensionsAnalysis an \code{\link{AFMImageFractalDimensionsAnalysis}} to store the results of the fractal analysis
#' @return a list of \code{\link{AFMImageFractalDimensionMethod}} objects with the calculated fractal dimensions and scales
#' @references  Gneiting2012, Tilmann Gneiting, Hana Sevcikova and Donald B. Percival 'Estimators of Fractal Dimension: Assessing the Roughness of Time Series and Spatial Data - Statistics in statistical Science, 2012, Vol. 27, No. 2, 247-277'
#' @author M.Beauvais
#' @rdname AFMFractalDimensionAnalyser-getFractalDimensions
#' @export
#' @seealso \code{\link{fractaldim}}
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' print(getFractalDimensions(AFMImageOfAluminiumInterface))
#
getFractalDimensions<-function(AFMImage, AFMImageFractalDimensionsAnalysis) {
  if (missing(AFMImageFractalDimensionsAnalysis)) {
    AFMImageFractalDimensionsAnalysis<-NULL
  }
  graphicalUpdate<-FALSE
  if (!is.null(AFMImageFractalDimensionsAnalysis)&&
      !is.null(AFMImageFractalDimensionsAnalysis@updateProgress)&&
      is.function(AFMImageFractalDimensionsAnalysis@updateProgress)&&
      !is.null(AFMImageFractalDimensionsAnalysis@updateProgress())) {
    graphicalUpdate<-TRUE
    
  }
  if (graphicalUpdate) {
    AFMImageFractalDimensionsAnalysis@updateProgress(message="1/2 - Calculating table", value=0)
  }
  
  totalLength <- 8
  counter<-0
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  rf2d <- matrix(AFMImage@data$h, nrow=AFMImage@samplesperline)
  fullfilename<-AFMImage@fullfilename
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  fd2d_transectvar <- fd.estim.transect.var(rf2d, p.index = 1, direction='hvd+d-')
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  fd2d_transectincr1 <- fd.estim.transect.incr1(rf2d, p.index = 1, direction='hvd+d-')
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  fd2d_isotropic <- fd.estim.isotropic(rf2d, p.index = 1, direction='hvd+d-')
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  fd2d_squareincr <- fd.estim.squareincr(rf2d, p.index = 1)
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  fd2d_filter1 <- fd.estim.filter1(rf2d, p.index = 1)
  
  if (graphicalUpdate) {
    counter<-counter+1
    value<-counter / totalLength
    text <- paste0(round(counter, 2),"/",totalLength)
    AFMImageFractalDimensionsAnalysis@updateProgress(value= 0, detail = text)
  }
  res=c(new("AFMImageFractalDimensionMethod", fd_method = "isotropic", fd = fd2d_isotropic$fd , fd_scale = fd2d_isotropic$scale))
  res=c(res, new("AFMImageFractalDimensionMethod", fd_method = "transectvar", fd = fd2d_transectvar$fd, fd_scale = fd2d_transectvar$scale))
  res=c(res, new("AFMImageFractalDimensionMethod", fd_method = "transectincr1", fd = fd2d_transectincr1$fd, fd_scale = fd2d_transectincr1$scale))
  res=c(res, new("AFMImageFractalDimensionMethod", fd_method = "squareincr", fd = fd2d_squareincr$fd, fd_scale = fd2d_squareincr$scale))
  res=c(res, new("AFMImageFractalDimensionMethod", fd_method = "filter1", fd = fd2d_filter1$fd, fd_scale = fd2d_filter1$scale))
  
  return(res)
}

exportFractalDimImagesForReport<-function(AFMImage, reportDirectory) {
  warning("Possible inconsistency between fractaldim images and values")
  sampleName<-basename(AFMImage@fullfilename)
  rf2d <- matrix(AFMImage@data$h, nrow=AFMImage@samplesperline)
  
  png(getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "isotropic"))
  fd2d_isotropic <- fd.estim.isotropic(rf2d, p.index = 1, direction='hvd+d-', plot.loglog = TRUE, plot.allpoints = TRUE)
  dev.off()
  png(getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "squareincr"))
  fd2d_squareincr <- fd.estim.squareincr(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)
  dev.off()
  png(getFractalDimensionsPngFullfilename(reportDirectory, sampleName, "filter1"))
  fd2d_filter1 <- fd.estim.filter1(rf2d, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)
  dev.off()
}

getFractalDimensionsPngFullfilename<-function(reportDirectory, imagebasename, method) {
  return(paste(reportDirectory, paste(imagebasename, "-", method, "-fractaldim.png", sep=""), sep="/"))
}