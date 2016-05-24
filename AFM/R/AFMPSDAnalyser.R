require("fftwtools")
require("pracma")

require("data.table")

require("gstat")
require(sp)

require("stringr")

# normality tests
require(gridExtra)
require(ggplot2)

require(reshape2)

require(stats)

if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("r", "roughness","x","predict.gstat"))

#' @title AFM image Power Spectrum Density analysis class
#' 
#' @description \code{AFMImagePSDAnalysis} handles an \code{\link{AFMImage}} roughness against lenghscale analysis 
#'
#' @slot roughnessAgainstLengthscale a data.table to store the roughness against lengthscale data
#' @slot intersections a list to store the lengthscales values as the intersections between slopes and the sill in roughness against lenghscale graph
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
#' @author M.Beauvais
AFMImagePSDAnalysis<-setClass("AFMImagePSDAnalysis",
                              slots = c(
                                psd1d_breaks="numeric",
                                psd2d_truncHighLengthScale="logical",
                                psd2d_maxHighLengthScale="numeric",
                                psd1d="data.table",
                                psd2d="data.table",
                                roughnessAgainstLengthscale="data.table",
                                intersections="numeric",
                                updateProgress="function"),
                              validity = function(object) { 
                                return(TRUE)
                              }
)

#' Constructor method of AFMImagePSDAnalysis Class.
#' 
#' @param .Object an AFMImagePSDAnalysis object
#' @rdname AFMImagePSDAnalysis-class
#' @export
setMethod("initialize",
          "AFMImagePSDAnalysis",
          function(.Object) {
            .Object@psd1d_breaks<-32
            .Object@psd2d_truncHighLengthScale<-TRUE
            .Object@psd2d_maxHighLengthScale<-0
            .Object@psd1d<-data.table()
            .Object@psd2d<-data.table()
            .Object@roughnessAgainstLengthscale<-data.table()
            validObject(.Object) ## valide l'objet
            return(.Object)
          })

#' Wrapper function AFMImagePSDAnalysis
#'
#' @rdname AFMImagePSDAnalysis-class
#' @export
AFMImagePSDAnalysis <- function() {
  return(new("AFMImagePSDAnalysis"))
}


#' Method \code{psd1d_breaks} returns a number of breaks to calculate PSD1D from PSD2D
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("psd1d_breaks",function(object){standardGeneric("psd1d_breaks")})
setGeneric(name= "psd1d_breaks<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("psd1d_breaks<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases psd1d_breaks
#' @param object a \code{\link{AFMImagePSDAnalysis}}
setMethod("psd1d_breaks",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@psd1d_breaks)
          }
)
setReplaceMethod(f="psd1d_breaks",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "numeric"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@psd1d_breaks <- value
                   return(AFMImagePSDAnalysis)
                 })

#' Method \code{psd2d_maxHighLengthScale} returns the maximum lengthscale to be managed by PSD2D
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("psd2d_maxHighLengthScale",function(object){standardGeneric("psd2d_maxHighLengthScale")})
setGeneric(name= "psd2d_maxHighLengthScale<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("psd2d_maxHighLengthScale<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases psd2d_maxHighLengthScale
setMethod("psd2d_maxHighLengthScale",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@psd2d_maxHighLengthScale)
          }
)
setReplaceMethod(f="psd2d_maxHighLengthScale",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "numeric"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@psd2d_maxHighLengthScale <- value
                   return(AFMImagePSDAnalysis)
                 })

#' Method \code{psd2d_truncHighLengthScale} returns if the lengthscale of PSD2D should be truncated
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("psd2d_truncHighLengthScale",function(object){standardGeneric("psd2d_truncHighLengthScale")})
setGeneric(name= "psd2d_truncHighLengthScale<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("psd2d_truncHighLengthScale<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases psd2d_truncHighLengthScale
setMethod("psd2d_truncHighLengthScale",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@psd2d_truncHighLengthScale)
          }
)
setReplaceMethod(f="psd2d_truncHighLengthScale",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "logical"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@psd2d_truncHighLengthScale <- value
                   return(AFMImagePSDAnalysis)
                 })




#' Method \code{psd1d} returns a data.table of psd in 1D
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("psd1d",function(object){standardGeneric("psd1d")})
setGeneric(name= "psd1d<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("psd1d<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases psd1d
setMethod("psd1d",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@psd1d)
          }
)
setReplaceMethod(f="psd1d",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "data.table"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@psd1d <- value
                   return(AFMImagePSDAnalysis)
                 })


#' Method \code{psd2d} returns a data.table of psd in 1D
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("psd2d",function(object){standardGeneric("psd2d")})
setGeneric(name= "psd2d<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("psd2d<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases psd2d
setMethod("psd2d",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@psd2d)
          }
)
setReplaceMethod(f="psd2d",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "data.table"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@psd2d <- value
                   return(AFMImagePSDAnalysis)
                 })

#' Method \code{roughnessAgainstLengthscale} returns a data.table of roughnesses versus lengthscale
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("roughnessAgainstLengthscale",function(object){standardGeneric("roughnessAgainstLengthscale")})
setGeneric(name= "roughnessAgainstLengthscale<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("roughnessAgainstLengthscale<-"))
           })


#' @rdname AFMImagePSDAnalysis-class
#' @aliases roughnessAgainstLengthscale
setMethod("roughnessAgainstLengthscale",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@roughnessAgainstLengthscale)
          }
)
setReplaceMethod(f="roughnessAgainstLengthscale",
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "data.table"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@roughnessAgainstLengthscale <- value
                   return(AFMImagePSDAnalysis)
                 })

#' Method \code{intersections} returns a intersection numeric value
#' @name AFMImagePSDAnalysis-class
#' @rdname AFMImagePSDAnalysis-class
setGeneric("intersections",function(object){standardGeneric("intersections")})
setGeneric(name= "intersections<-", 
           def= function(AFMImagePSDAnalysis, value) {
             return(standardGeneric("intersections<-"))
           })

#' @rdname AFMImagePSDAnalysis-class
#' @aliases intersections
setMethod("intersections",signature=signature(object='AFMImagePSDAnalysis'),
          function(object) {
            return(object@intersections)
          }
)
setReplaceMethod(f="intersections", 
                 signature(AFMImagePSDAnalysis = "AFMImagePSDAnalysis", value = "numeric"),
                 definition= function(AFMImagePSDAnalysis, value) {
                   AFMImagePSDAnalysis@intersections <- value
                   return(AFMImagePSDAnalysis)
                 })


#' Shift the quadrants of the FFT 2D
#'
#' \code{shiftFFT2D} returns the FFT 2D matrix shifted to put zero frequencies in the middle.
#' 
#' @param fft2data  the FFT 2D of the AFM image
#' @return The shifted matrix
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' library(fftwtools)
#' 
#' data(AFMImageOfNormallyDistributedHeights)
#' AFMImage<-AFMImageOfNormallyDistributedHeights
#' nMheightsData= matrix(AFMImage@@data$h, nrow=AFMImage@@samplesperline)
#' shiftedFFT2D<-shiftFFT2D(fftw2d(nMheightsData))
shiftFFT2D<-function(fft2data) { 
  N=nrow(fft2data)
  M=ncol(fft2data)
  halfN=N/2
  halfM=M/2
  quadrant1=fft2data[1:halfN, seq(1,halfM)]
  quadrant2=fft2data[seq(halfN+1,N), seq(1,halfM)]
  quadrant3=fft2data[seq(halfN+1,N),seq(halfM+1,M)]
  quadrant4=fft2data[seq(1,halfN),seq(halfM+1,M)]
  return(rbind(cbind(quadrant3,quadrant2),cbind(quadrant4, quadrant1)))
}    


# zeroPadShiftedFFT2D<-function(shiftedFFT2Ddata){
#   N=nrow(fft2data)
#   M=ncol(fft2data)
#   
#   r = 2^(ceil(log2(x)))
# }

#' Calculate the shifted PSD matrix
#'
#' \code{shiftedPSDuv} returns the Power Spectral Density matrix in the frequency space from shifted FFT 2D
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @return (1/NM^2) * abs(shiftedFFT2Ddata)^2) with N the number of lines of the sample and M the number of samples per line of the sample
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' library(ggplot2)
#' 
#' data(AFMImageOfRegularPeaks)
#' AFMImage<-AFMImageOfRegularPeaks
#' nMheightsData= matrix(AFMImage@@data$h, nrow=AFMImage@@samplesperline)
#' shiftedPSDuv<-shiftedPSDuv(AFMImage)

#' a=AFMImage@@scansize
#' b=AFMImage@@scansize
#' 
#' M=AFMImage@@sampsline
#' N=AFMImage@@lines
#' NM=N*M # pixels^2
#' MN = M*N 
#' A=a*b
#' ab=a*b
#' 
#' dx=a/M
#' dy=b/N
#' 
#' um = seq( (1-(M+1)/2)/(M*dx), (M-(M+1)/2)/(M*dx), by=1/(M*dx))
#' vn = seq( (1-(N+1)/2)/(N*dy), (N-(N+1)/2)/(N*dy), by=1/(N*dy))
#' x = rep(um, times = AFMImage@@lines)
#' y = rep(vn, each = AFMImage@@sampsline)
#' z = as.vector(shiftedPSDuv)
#' 
#' data<-data.frame(x=x, y=y, z=z)
#' 
#' p5 <- qplot(x, y, data=data, colour=log10(z))
#' p5 <- p5 + scale_colour_gradientn(colours = rainbow(7))
#' p5 <- p5 + ylab("v")
#' p5 <- p5 + xlab("u")
#' title<-paste("shifted PSD of", basename(AFMImage@@fullfilename))
#' p5 <- p5 + ggtitle(title)
#' # Hide all the horizontal gridlines
#' p5 <- p5 + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
#' # Hide all the vertical gridlines
#' p5 <- p5 + theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
#' p5 <- p5 + theme(panel.background = element_rect(fill = 'white', colour = 'black'))
#' p5
shiftedPSDuv<-function(AFMImage) {
  nMheighData= matrix(AFMImage@data$h, nrow=AFMImage@samplesperline)
  shiftedFFT2Ddata = shiftFFT2D(fftw2d(nMheighData))
  
  N=nrow(shiftedFFT2Ddata)
  M=ncol(shiftedFFT2Ddata)
  NM=N*M
  return((1/NM^2) * abs(shiftedFFT2Ddata)^2)
}

#' Calculate the 2D Power Spectral Density
#' 
#' PSD2DAgainstFrequency returns a data table of PSD 2D values against spatial frequencies
#'
#' @param AFMImage an \code{AFMImage} to be analysed
#' @param AFMImagePSDAnalysis an \code{AFMImagePSDAnalysis} to store PSD analysis results
#' @return \code{PSD2DAgainstFrequency} returns a data table of frequencies and PSD values
#' \itemize{
#'   \item freq: the considered frequency
#'   \item PSD: the considered PSD value
#'   \item type: PSD-2D
#'   \item fullfilename: directory and filename on the disk
#' }
#' @references  Sidick2009, Erkin Sidick "Power Spectral Density Specification and Analysis of Large Optical Surfaces", 2009, "Modeling Aspects in Optical Metrology II, Proc. of SPIE Vol. 7390 73900L-1"
#' @name PSD2DAgainstFrequency
#' @rdname PSD2DAgainstFrequency-methods
#' @exportMethod PSD2DAgainstFrequency
#' @examples
#' \dontrun{
#' library(AFM)
#' library(ggplot2)
#' library(plyr)
#' 
#' # Calculate Power Spectrum Density in 2D against frequency
#' data("AFMImageOfNormallyDistributedHeights")
#' oneAFMImage<-AFMImageOfNormallyDistributedHeights
#' psd2d<-PSD2DAgainstFrequency(oneAFMImage)
#' p <- ggplot(data=psd2d)
#' p <- p + geom_point(aes(freq, PSD, color=type),subset = .(type %in% c("PSD-2D")))
#' p <- p + geom_line(aes(freq, PSD, color=type),subset = .(type %in% c("PSD-1D")),size=1.1)
#' p <- p + scale_x_log10()
#' p <- p + scale_y_log10()
#' p <- p + ylab("PSD (nm^4)")
#' p <- p + xlab("Frequency (nm^-1)")
#' p <- p + ggtitle(basename(oneAFMImage@@fullfilename))
#' p
#' }
setGeneric(name= "PSD2DAgainstFrequency", 
           def= function(AFMImage, AFMImagePSDAnalysis) {
             return(standardGeneric("PSD2DAgainstFrequency"))
           })

#' @rdname PSD2DAgainstFrequency-methods
#' @aliases PSD2DAgainstFrequency,AFMImage-method
setMethod(f="PSD2DAgainstFrequency", "AFMImage",
          definition= function(AFMImage, AFMImagePSDAnalysis) {
            NyquistFq<-getNyquistSpatialFrequency(AFMImage)
            
            a=AFMImage@hscansize
            b=AFMImage@vscansize
            
            M=AFMImage@samplesperline
            N=AFMImage@lines
            NM=N*M # pixels^2
            MN = M*N 
            A=a*b
            ab=a*b
            
            dx=a/M
            dy=b/N
            
            shiftedPSDuv<-shiftedPSDuv(AFMImage)
            
            um = seq( (1-(M+1)/2)/(M*dx), (M-(M+1)/2)/(M*dx), by=1/(M*dx))
            vn = seq( (1-(N+1)/2)/(N*dy), (N-(N+1)/2)/(N*dy), by=1/(N*dy))
            K=meshgrid(um,vn)
            K$Z<-sqrt(K$X^2+K$Y^2)
            
            aAggregatedPSDValuesForEachFreq=data.frame(freq=sort(unique(as.vector(K$Z))))
            totalLength <- length(aAggregatedPSDValuesForEachFreq$freq)
            frequencies<-c()
            sumedPSD<-c()
            counter<-0
            for(freq in aAggregatedPSDValuesForEachFreq$freq) {
              if (freq > NyquistFq) break;
              
              if (!is.null(AFMImagePSDAnalysis@updateProgress)&&
                  (is.function(AFMImagePSDAnalysis@updateProgress)&&
                   (!is.null(AFMImagePSDAnalysis@updateProgress())))) {
                counter<-counter+1
                if (counter/100==floor(counter/100)) {
                  value<-counter / totalLength
                  text <- paste0("freq:", round(freq, 2)," ", round(counter, 2),"/",totalLength)
                  AFMImagePSDAnalysis@updateProgress(value= value, detail = text)
                } 
              }
              
              inds <- arrayInd(which(K$Z == freq), dim(K$Z))
              allPSDSum<-0
              allPSDSum<-sum(shiftedPSDuv[inds[,1:2]])
              sumedPSD = c(sumedPSD, allPSDSum)
              frequencies=c(frequencies, freq)
            }
            return(data.table(freq = frequencies, PSD = sumedPSD, type="PSD-2D", name=AFMImage@fullfilename))
          }
)

#' Calculate the 1D Power Spectral Density; returns a data table of PSD 1D and PSD 2D values
#' against spatial frequencies.\cr As mentionned in Sidick2009, this function calculates the 
#' PSD against spatial frequencies in 1D from \code{\link{PSD2DAgainstFrequency}} by using
#' breaks in the log space to sum PSD 2D and frequency values.
#' 
#' @param AFMImage an \code{AFMImage} to be analysed
#' @param AFMImagePSDAnalysis n \code{AFMImagePSDAnalysis} to store the setup and results of PSD analysis
#' @return \code{PSD1DAgainstFrequency} returns a data table of frequencies and PSD values
#' \itemize{
#'   \item freq: the considered frequency
#'   \item PSD: the considered PSD value
#'   \item type: PSD-1D
#'   \item fullfilename: directory and filename on the disk
#' }
#' @name PSD1DAgainstFrequency
#' @rdname PSD1DAgainstFrequency-methods
#' @exportMethod PSD1DAgainstFrequency
#' @examples
#' \dontrun{
#' library(AFM)
#' library(ggplot2)
#' library(plyr)
#' library(scales)

#' data("AFMImageOfNormallyDistributedHeights")
#'  newAFMImage<-AFMImageOfNormallyDistributedHeights
#' newAFMImage@fullfilename<-"C:/Users/one/AFMImageOfNormallyDistributedHeights.txt"
#' psdAnalysis<-AFMImagePSDAnalysis()
#' # Create a closure to update progress
#' psdAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
#'   if (exists("progressPSD")){
#'    if (!is.null(message)) {
#'      progressPSD$set(message = message, value = 0)
#'    }else{
#'      progressPSD$set(value = value, detail = detail)
#'    }
#'   }
#' }
#' psdAnalysis@psd1d_breaks<-2^3
#' psdAnalysis@psd2d_truncHighLengthScale<-TRUE
#' psdAnalysis<-performAllPSDCalculation(AFMImagePSDAnalysis= psdAnalysis, AFMImage= newAFMImage)
#' datap<-psdAnalysis@psd1d
#' p <- ggplot(data=datap)
#' p <- p + geom_point(aes(freq, PSD, color=type),data=datap[datap$type %in% c("PSD-2D")])
#' p <- p + geom_line(aes(freq, PSD, color=type),data=datap[datap$type %in% c("PSD-1D")],size=1.1)
#' p <- p + scale_x_log10()
#' p <- p + scale_y_log10()
#' p <- p + ylab("PSD (nm^4)")
#' p <- p + xlab("Frequency (nm^-1)")
#' p
#' }
setGeneric(name= "PSD1DAgainstFrequency", 
           def= function(AFMImage,AFMImagePSDAnalysis) {
             return(standardGeneric("PSD1DAgainstFrequency"))
           })

#' @rdname PSD1DAgainstFrequency-methods
#' @aliases PSD1DAgainstFrequency,AFMImage-method
setMethod(f="PSD1DAgainstFrequency", "AFMImage",
          definition= function(AFMImage, AFMImagePSDAnalysis) {
            AFMImagePSDAnalysis@psd2d<-PSD2DAgainstFrequency(AFMImage, AFMImagePSDAnalysis)
            
            breaks=AFMImagePSDAnalysis@psd1d_breaks
            psd2dDT=AFMImagePSDAnalysis@psd2d
            
            # step 3, cut in the log space
            Q <- breaks
            maxRhoL<- max(psd2dDT$freq)
            maxRhoL
            
            psd2dDT$logcuts<-cut(log10(psd2dDT$freq),breaks = Q)
            
            meanFreq<-c()
            meanPSD<-c()
            totalLength<-length(unique(as.vector(psd2dDT$logcuts)))
            counter<-0
            if (!is.null(AFMImagePSDAnalysis@updateProgress)&&
                (is.function(AFMImagePSDAnalysis@updateProgress)&&
                 (!is.null(AFMImagePSDAnalysis@updateProgress())))) {
              text <- paste0("starting ", totalLength, " calculations")
              AFMImagePSDAnalysis@updateProgress(value= 0, detail = text)
            } 
            
            
            for(freq in sort(unique(as.vector(psd2dDT$logcuts)))) {
              inds <- arrayInd(which(psd2dDT$logcuts == freq), dim(psd2dDT))
              allFreqSum<-0
              allPSDSum<-0
              allFreqSum<-mean(psd2dDT$freq[inds[,1]])
              allPSDSum<-mean(psd2dDT$PSD[inds[,1]])
              meanFreq = c(meanFreq, allFreqSum)
              meanPSD = c(meanPSD, allPSDSum)
              
              if (!is.null(AFMImagePSDAnalysis@updateProgress)&&
                  (is.function(AFMImagePSDAnalysis@updateProgress)&&
                   (!is.null(AFMImagePSDAnalysis@updateProgress())))) {
                counter<-counter+1
                if (counter/100==floor(counter/100)) {
                  value<- counter / totalLength
                  text <- paste0("freq:", round(freq, 2)," ", round(counter, 2),"/",totalLength)
                  AFMImagePSDAnalysis@updateProgress(value= value, detail = text)
                } 
              }
            }
            return(rbind(data.table(freq = meanFreq, PSD = meanPSD, type="PSD-1D", name=AFMImage@fullfilename), data.table(freq = psd2dDT$freq, PSD = psd2dDT$PSD, type="PSD-2D", name=psd2dDT$name)))
          })

#' Calculate the roughness of the sample against length scale
#'
#' The calculation of the roughness against lengthscale is performed throught a FFT 2D calculation, PSD 2D calculation and a meshgrid of frequencies.
#' \code{RoughnessByLengthScale} returns a data.table of roughnesses against length scales
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImagePSDAnalysis n \code{AFMImagePSDAnalysis} to store the setup and results of PSD analysis
#' 
#' @return a data table of lenght scale (r) and roughness values (roughness)
#' \itemize{ 
#' \item {roughness: roughnesses} 
#' \item {r: length scales}
#' \item {filename: fullfilename slot of the AFMImage} 
#' }
#' @name RoughnessByLengthScale
#' @rdname RoughnessByLengthScale-methods
#' @exportMethod RoughnessByLengthScale
#' @author M.Beauvais
#' @examples
#' library(AFM)
#' library(ggplot2)
#' 
#' data("AFMImageOfNormallyDistributedHeights")
#' oneAFMImage<-AFMImageOfNormallyDistributedHeights
#' AFMImagePSDAnalysis<-AFMImagePSDAnalysis()
#' data<-RoughnessByLengthScale(oneAFMImage, AFMImagePSDAnalysis)
#' r<-roughness<-filename<-NULL
#' p1 <- ggplot(data, aes(x=r, y=roughness, colour= basename(filename)))
#' p1 <- p1 + geom_point()
#' p1 <- p1 + geom_line()
#' p1 <- p1 + ylab("roughness (nm)")
#' p1 <- p1 + xlab("lengthscale (nm)")
#' p1
setGeneric(name= "RoughnessByLengthScale", 
           def= function(AFMImage, AFMImagePSDAnalysis) {
             return(standardGeneric("RoughnessByLengthScale"))
           })

#' @rdname RoughnessByLengthScale-methods
#' @aliases RoughnessByLengthScale,AFMImage-method
setMethod(f="RoughnessByLengthScale", "AFMImage",
          definition= function(AFMImage, AFMImagePSDAnalysis) {
            # calculate roughness depending on frequency
            AFMImagePSDAnalysis@psd2d<-PSD2DAgainstFrequency(AFMImage, AFMImagePSDAnalysis)
            
            
            truncHighLengthScale = AFMImagePSDAnalysis@psd2d_truncHighLengthScale
            maxHighLengthScale = AFMImagePSDAnalysis@psd2d_maxHighLengthScale
            AggregatedPSDValuesForEachFreq = AFMImagePSDAnalysis@psd2d
            
            minFrequency<-1/min(AFMImage@hscansize, AFMImage@vscansize)
            indexfmin<-tail(which(AggregatedPSDValuesForEachFreq$freq < minFrequency), n=1)
            
            if (missing(truncHighLengthScale)||truncHighLengthScale==FALSE) {
              if(!missing(maxHighLengthScale)){
                truncHighLengthScale <- FALSE
                if (maxHighLengthScale<(1/minFrequency)) {
                  indexfmin<-which(AggregatedPSDValuesForEachFreq$freq > (1/maxHighLengthScale))[1]-1
                }
              }
            }
            
            nyquistSF <- getNyquistSpatialFrequency(AFMImage)
            indexfmax<-which(AggregatedPSDValuesForEachFreq$freq > nyquistSF)[1]-1
            
            #if (!isTRUE(truncHighLengthScale)||is.na(indexfmin)) indexfmin<-0
            if (is.na(indexfmin)) indexfmin<-0
            if (is.na(indexfmax)) indexfmax<-length(AggregatedPSDValuesForEachFreq$freq)
            
            r<-c()
            roughnesses=c()
            totalLength<-indexfmax
            counter<-0
            for (i in seq(1,indexfmax)){
              if (i>indexfmin) {
                tryingPSDSum<-sum(AggregatedPSDValuesForEachFreq$PSD[i:indexfmax])
                roughnesses=c(roughnesses, sqrt(tryingPSDSum))
                r=c(r, 1/AggregatedPSDValuesForEachFreq$freq[i])
                
                if (!is.null(AFMImagePSDAnalysis@updateProgress)&&
                    is.function(AFMImagePSDAnalysis@updateProgress)&&
                    !is.null(AFMImagePSDAnalysis@updateProgress())) {
                  counter<-counter+1
                  if (counter/100==floor(counter/100)) {
                    value<-counter / totalLength
                    text <- paste0(round(counter, 2),"/",totalLength)
                    AFMImagePSDAnalysis@updateProgress(value= value, detail = text)
                  } 
                }
              }
            }
            return(data.table(filename=rep(AFMImage@fullfilename, length(AggregatedPSDValuesForEachFreq$freq)-indexfmin), r= r, roughness= roughnesses))
          })


#' Get the Nyquist spatial frequency
#'
#' Get the Nyquist spatial frequency of an \code{\link{AFMImage}} calculated as following:\cr
#' 0.5 multiplied by the minimum between the horizontal scansize divided by the number of samples per line and the vertical scansize divided by the number of lines
#' 
#' \code{getNyquistSpatialFrequency} returns the Nyquist spatial frequency as a numeric
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @return the Nyquist spatial frequency of the \code{\link{AFMImage}}
#' @name getNyquistSpatialFrequency
#' @rdname getNyquistSpatialFrequency-methods
#' @exportMethod getNyquistSpatialFrequency
#' @author M.Beauvais
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfNormallyDistributedHeights)
#' NyquistSpatialFrequency<-getNyquistSpatialFrequency(AFMImageOfNormallyDistributedHeights)
#' print(NyquistSpatialFrequency)
#' 
setGeneric(name= "getNyquistSpatialFrequency", 
           def= function(AFMImage) {
             return(standardGeneric("getNyquistSpatialFrequency"))
           })

#' @rdname getNyquistSpatialFrequency-methods
#' @aliases getNyquistSpatialFrequency,AFMImage-method
setMethod(f="getNyquistSpatialFrequency", "AFMImage",
          definition= function(AFMImage) {
            M=AFMImage@samplesperline
            N=AFMImage@lines
            a=AFMImage@hscansize
            b=AFMImage@vscansize
            dx=a/M
            dy=b/N
            
            #old   return(min(abs((1-(M+1)/2)/(M*dx)), abs((1-(N+1)/2)/(N*dy))))
            return(min(1/(2*dx),1/(2*dy)))
          })

#' Get a zero padded AFMImage
#' 
#' Get a zero padded \code{\link{AFMImage}} useful in Power Spectral Density analysis. The original \code{\link{AFMImage}} is padded with zero in order to get a larger square AFMImage which size is a power of 2.
#'
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @return a zero-padded \code{\link{AFMImage}} with a fullfilename equals to the original fullfilename pasted with padded-to-"ScanSize".txt
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfNormallyDistributedHeights)
#' paddedAFMImage<-getPaddedAFMImage(AFMImageOfNormallyDistributedHeights)
#' displayIn3D(AFMImage= paddedAFMImage, width= 1024)
getPaddedAFMImage<-function(AFMImage) {
  paddedAFMImageMatrix<-matrix(AFMImage@data$h, nrow=AFMImage@samplesperline, ncol=AFMImage@lines,byrow = TRUE)
  N=nrow(paddedAFMImageMatrix)
  print(N)
  M=ncol(paddedAFMImageMatrix)
  print(M)
  
  rn = 2^(ceil(log2(N)))
  paddedN <- ifelse(rn==N, 2^(ceil(log2(N+1))), rn)
  rm = 2^(ceil(log2(M)))
  paddedM <- ifelse(rm==M, 2^(ceil(log2(M+1))), rm)
  
  
  addingN=paddedN/4
  addingM=paddedM/4
  
  A<-matrix( rep(0,addingM*N), nrow=N,ncol=addingM,byrow = TRUE)
  B<-matrix( rep(0,addingM*paddedN), nrow=addingM,ncol=paddedN,byrow = TRUE)
  paddedAFMImageMatrix<-cbind(A, paddedAFMImageMatrix, A)
  paddedAFMImageMatrix<-rbind(B, paddedAFMImageMatrix, B)
  
  Lines<-paddedN;
  Samplesperline<-paddedM;
  ScanSize<-AFMImage@hscansize*paddedM/M
  # not tested
  hscanSize<-AFMImage@hscansize*paddedM/M
  vscanSize<-AFMImage@vscansize*paddedN/N
  ScanSize<-max(hscanSize, vscanSize)
  
  
  scanSizeFromZero<-ScanSize-1
  scanby<-ScanSize/Samplesperline
  endScan<-ScanSize*(1-1/Samplesperline)
  nM<-as.vector(t(paddedAFMImageMatrix))
  
  AFMImage(data = data.table(x = rep(seq(0,endScan, by= scanby), times = Lines),
                             y = rep(seq(0,endScan, by= scanby), each = Samplesperline), 
                             h = nM), 
           samplesperline = Samplesperline, 
           lines = Lines,
           hscansize = hscanSize,
           vscansize = vscanSize,
           scansize = ScanSize, 
           fullfilename = paste(AFMImage@fullfilename, "padded-to-",ScanSize,".txt",sep=""))
}

#' Perform all the calculation for PSD exploitation
#' 
#' \code{\link{performAllPSDCalculation}} perform all the calculation for PSD exploitation
#' @param AFMImagePSDAnalysis an \code{\link{AFMImagePSDAnalysis}} to manage and store the results of PSD analysis
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' data(AFMImageOfNormallyDistributedHeights)
#' 
#' newAFMImage<-AFMImageOfNormallyDistributedHeights
#' newAFMImage@fullfilename<-"C:/Users/one/AFMImageOfNormallyDistributedHeights.txt"
#' psdAnalysis<-AFMImagePSDAnalysis()
#' # Create a closure to update progress
#' psdAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
#'   if (exists("progressPSD")){
#'     if (!is.null(message)) {
#'       progressPSD$set(message = message, value = 0)
#'     }else{
#'       progressPSD$set(value = value, detail = detail)
#'     }
#'   }
#' }
#' psdAnalysis@psd1d_breaks<-2^3
#' psdAnalysis@psd2d_truncHighLengthScale<-TRUE
#' psdAnalysis<-performAllPSDCalculation(AFMImagePSDAnalysis= psdAnalysis, AFMImage= newAFMImage)
#' print("done psdAnalysis")
#' }
performAllPSDCalculation<-function(AFMImagePSDAnalysis, AFMImage) {
  if (is.function(AFMImagePSDAnalysis@updateProgress)) {
    AFMImagePSDAnalysis@updateProgress(message="1/3 - Calculating PSD2D", value=0)
  }
  AFMImagePSDAnalysis@psd2d<-PSD2DAgainstFrequency(AFMImage, AFMImagePSDAnalysis)
  
  if (is.function(AFMImagePSDAnalysis@updateProgress)) {
    AFMImagePSDAnalysis@updateProgress(message="2/3 Calculating PSD1D", value=0)
  }
  AFMImagePSDAnalysis@psd1d<-PSD1DAgainstFrequency(AFMImage, AFMImagePSDAnalysis)
  
  if (is.function(AFMImagePSDAnalysis@updateProgress)) {
    AFMImagePSDAnalysis@updateProgress(message="3/3 Calculating Roughness", value=0)
  }
  AFMImagePSDAnalysis@roughnessAgainstLengthscale<-RoughnessByLengthScale(AFMImage, AFMImagePSDAnalysis)
  
  return(AFMImagePSDAnalysis)
}

saveOnDiskIntersectionForRoughnessAgainstLengthscale<-function(AFMImageAnalyser, exportDirectory){
  sampleName<-basename(AFMImageAnalyser@fullfilename)
  
  data<-AFMImageAnalyser@psdAnalysis@roughnessAgainstLengthscale
  data$r<-as.numeric(data$r)
  
  aval<-max(data$r)
  index<-which(data$r<= aval)[1]
  
  lengthData<-length(data$r)-index
  ndataw<-tail(data,n= lengthData)
  ndataw$sample<-basename(ndataw$filename)
  
  #find x1 x2 that minimizes Xinter
  min=1
  max=2
  
  point<- data[data$r %in% min(data$r)]
  otherpoints<-as.numeric(point$V1)-min
  point1<-data[as.numeric(data$V1) %in% (otherpoints)]
  otherpoints<-as.numeric(point$V1)-max
  point2<-data[as.numeric(data$V1) %in% (otherpoints)]
  
  point1<-data[min]
  point2<-data[max]
  
  
  origintangeantePoints = data.table(x=c(point1$r, point2$r), y=c(point1$roughness, point2$roughness))
  aorigin<-0
  borigin <- point1$roughness - aorigin * point2$r
  
  
  
  x1x2=AFMImageAnalyser@psdAnalysis@intersections[c(2,3,5,6)]
  
  
  for(i in c(0,2)) {
    x1<-x1x2[i+1]
    x2<-x1x2[i+2]
    #x1<-244
    #x2<-44
    point1<-data[x1]
    point2<-data[x2]
    x=data[seq(x1,x2)]$r
    y=data[seq(x1,x2)]$roughness
    res <- lm(y~x)
    coefficients(res)
    b<-unname(res$coefficients[1])
    a<-unname(res$coefficients[2])
    
    tangeantePoints = data.table(x=c(point1$r, point2$r), y=c(point1$roughness, point2$roughness))
    
    xinter <- (b-borigin)/(aorigin-a)
    title<-paste("Xinter= ", xinter," -plateau= ", borigin)
    
    roughness<-r<-NULL
    p1 <- ggplot(ndataw, aes(x=r, y=roughness, colour= basename(sample)))
    p1 <- p1 + geom_point()
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_abline(intercept = b, slope = a)
    p1 <- p1 +geom_point(data=tangeantePoints, aes(x=x, y=y), color="blue") 
    p1 <- p1 + geom_abline(intercept = borigin, slope = aorigin)
    p1 <- p1 +geom_point(data=origintangeantePoints, aes(x=x, y=y), color="blue") 
    p1 <- p1 + ylab("roughness (nm)")
    p1 <- p1 + xlab("lengthscale (nm)")
    p1 <- p1 + guides(colour=FALSE)
    p1 <- p1 + ggtitle(title)
    
    exportpng2FullFilename=getRoughnessAgainstLengthscaleIntersection(exportDirectory, paste( sampleName, i, sep="-"))
    print(paste("saving", basename(exportpng2FullFilename)))
    png(filename=exportpng2FullFilename, units = "px", width=800, height=800)
    print(p1)
    dev.off()
  }
}


getAutoIntersectionForRoughnessAgainstLengthscale<-function(psdAnalysis, AFMImage, second_slope=FALSE){
  sampleName<-basename(AFMImage@fullfilename)
  exportDirectory<-paste(dirname(AFMImage@fullfilename), "outputs", sep="/")
  
  data<-psdAnalysis@roughnessAgainstLengthscale
  data$r<-as.numeric(data$r)
  
  aval<-max(data$r)
  index<-which(data$r<= aval)[1]
  
  lengthData<-length(data$r)-index
  #print(paste("lengthData=",lengthData))
  #lengthData
  ndataw<-tail(data,n= lengthData)
  ndataw$sample<-basename(ndataw$filename)
  #print(paste("length(ndataw)=",length(ndataw)))
  
  
  
  minimumR <- function(data, space, x, y) {
    lengthData<-nrow(data)
    
    aorigin<-0
    borigin <- data[1]$roughness
    #print(borigin)
    
    finalres2=c()
    finalres = c(Inf,0,0)
    for (i in seq(1, length(x))) {
      x1=x[i]
      for (j in seq(1, length(y))) {
        if (abs(j-i)>space) {
          x2=y[j]
          
          if ((x1<1)||(x2<1)||(x1>lengthData)||(x2>lengthData)||(x1==x2)) {
            inter <- data[1]$r
          } else{
            if (x1<x2) {
              myx=data[seq(x1,x2)]$r
              myy=data[seq(x1,x2)]$roughness
            }
            if (x1>x2) {
              myx=data[seq(x2,x1)]$r
              myy=data[seq(x2,x1)]$roughness
            }    
            
            res <- lm(myy~myx)
            b<-unname(res$coefficients[1])
            a<-unname(res$coefficients[2])
            inter <- (borigin-b)/a
            if (inter<finalres[1]) {
              finalres=c(inter, x1, x2)
              #print(finalres)
              
            }
          }
        }
        #print(paste(x1, x2))
      }
      #finalres2=c(finalres2, inter)
    }
    return(finalres)
  }
  
  lengthData<-nrow(data)
  
  if (second_slope==FALSE) {
    aby<-floor(lengthData/10)
    x <- seq(1,lengthData,by=aby)
    z <- minimumR(data, space= 1, x,x)
  } else {
    newMax = lengthData-floor(lengthData*90/100)
    #print(newMax)
    aby<-ceil(newMax/150)
    #print(aby)
    if (newMax<100) space = 1
    else space = 100
    #space
    x <- seq(1,newMax,by=aby)
    z <- minimumR(data, space= space, x,x)
  }
  return(z)
}

getRoughnessAgainstLengthscale<-function(exportDirectory, sampleName) {
  exportCsvFilename<-paste(sampleName,"-roughness-against-lengthscale.png", sep="")
  exportCsvFullFilename<-paste(exportDirectory, exportCsvFilename, sep="/")
  return(exportCsvFullFilename)
}

getRoughnessAgainstLengthscale10nm<-function(exportDirectory, sampleName) {
  exportCsvFilename<-paste(sampleName,"-roughness-against-lengthscale-10nm.png", sep="")
  exportCsvFullFilename<-paste(exportDirectory, exportCsvFilename, sep="/")
  return(exportCsvFullFilename)
}

getRoughnessAgainstLengthscaleIntersection<-function(exportDirectory, sampleName) {
  exportpngFilename<-paste(sampleName,"-roughness-against-lengthscale-intersection.png", sep="")
  exportpngFullFilename<-paste(exportDirectory, exportpngFilename, sep="/")
  return(exportpngFullFilename)
}
