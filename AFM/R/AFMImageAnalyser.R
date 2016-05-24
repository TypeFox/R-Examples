require("fftwtools")
require("pracma")

require("data.table")

require("gstat")
require(sp)


require("stringr")

# normality tests
require(gridExtra)
require(moments)    
require(ggplot2)

require(reshape2)

# for reporting
require(png)
require(grid)

if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("h", "..density.."))

#' AFM image analyser class
#' 
#' A S4 class to handle the analysis of one AFM Image. 
#'
#' @slot AFMImage \code{\link{AFMImage}} to be analysed
#' @slot variogramAnalysis \code{\link{AFMImageVariogramAnalysis}}
#' @slot psdAnalysis \code{\link{AFMImagePSDAnalysis}}
#' @slot fdAnalysis \code{\link{AFMImageFractalDimensionsAnalysis}}
#' @slot networksAnalysis \code{\link{AFMImageNetworksAnalysis}}
#' @slot mean  the mean of heights of the \code{\link{AFMImage}}
#' @slot variance the variance of heights of the \code{\link{AFMImage}}
#' @slot TotalRrms the total Root Mean Square Roughness of the \code{\link{AFMImage}} calculated from variance
#' @slot Ra mean roughness or mean of absolute values of heights
#' @slot fullfilename to be removed ?
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImageAnalyser-class
#' @rdname AFMImageAnalyser-class
#' @exportClass AFMImageAnalyser 
#' @author M.Beauvais
#' @include AFMVariogramAnalyser.R AFMPSDAnalyser.R AFMFractalDimensionAnalyser.R AFM3DPrinter.R AFMNetworksAnalyser.R AFM3DPrinter.R
#'
AFMImageAnalyser<-setClass("AFMImageAnalyser",
                           slots = c(
                             versions="data.table",
                             AFMImage="AFMImage",
                             variogramAnalysis="AFMImageVariogramAnalysis", 
                             psdAnalysis="AFMImagePSDAnalysis",
                             fdAnalysis="AFMImageFractalDimensionsAnalysis",
                             networksAnalysis="AFMImageNetworksAnalysis",
                             threeDimensionAnalysis="AFMImage3DModelAnalysis",
                             mean="numeric", 
                             variance="numeric", 
                             TotalRrms="numeric", 
                             Ra="numeric", 
                             fullfilename="character",
                             updateProgress="function"))

#' Constructor method of AFMImageAnalyser Class.
#'
#' @param .Object an AFMImageAnalyser object
#' @param AFMImage an \code{AFMImage}
#' @param variogramAnalysis \code{\link{AFMImageVariogramAnalysis}}
#' @param psdAnalysis \code{\link{AFMImagePSDAnalysis}}
#' @param fdAnalysis \code{\link{AFMImageFractalDimensionsAnalysis}}
#' @param networksAnalysis \code{\link{AFMImageNetworksAnalysis}}
#' @param threeDimensionAnalysis \code{\link{AFMImage3DModelAnalysis}}
#' @param mean  the mean of heights of the \code{\link{AFMImage}}
#' @param variance the variance of heights of the \code{\link{AFMImage}}
#' @param TotalRrms the total Root Mean Square Roughness of the \code{\link{AFMImage}} calculated from variance
#' @param Ra mean roughness or mean of absolute values of heights
#' @param fullfilename to be removed?
#' @rdname AFMImageAnalyser-class-initialize
#' @export
setMethod("initialize", "AFMImageAnalyser", function(.Object,
                                                     AFMImage,
                                                     variogramAnalysis, 
                                                     psdAnalysis,
                                                     fdAnalysis,
                                                     networksAnalysis,
                                                     threeDimensionAnalysis,
                                                     mean, 
                                                     variance, 
                                                     TotalRrms, 
                                                     Ra, 
                                                     fullfilename)  
{
  if (!missing(AFMImage)) .Object@AFMImage<-AFMImage
  if (!missing(variogramAnalysis)) .Object@variogramAnalysis<-variogramAnalysis
  if (!missing(psdAnalysis)) .Object@psdAnalysis<-psdAnalysis
  if (!missing(fdAnalysis)) .Object@fdAnalysis<-fdAnalysis
  if (!missing(networksAnalysis)) .Object@networksAnalysis<-networksAnalysis
  if (!missing(threeDimensionAnalysis)) .Object@threeDimensionAnalysis<-threeDimensionAnalysis
  if (!missing(mean)) .Object@mean<-mean
  if (!missing(variance)) .Object@variance<-variance
  if (!missing(TotalRrms)) .Object@TotalRrms<- TotalRrms
  if (!missing(Ra)) .Object@Ra<-Ra
  .Object@fullfilename<-fullfilename
  .Object@versions<-getLibrariesVersions()
  validObject(.Object)      
  return(.Object)
})

#' Wrapper function AFMImageAnalyser
#'
#' @param AFMImage an \code{AFMImage}
#' @rdname AFMImageAnalyser-class
#' @export
AFMImageAnalyser <- function(AFMImage) {
  return(new("AFMImageAnalyser", AFMImage= AFMImage, fullfilename= AFMImage@fullfilename))
}

#' Analyse an AFMImage
#' 
#' A function to wrap all the analysis of an \code{\link{AFMImage}}
#' \itemize{
#'   \item variogram analysis  including evaluation of basic variogram models with sill and range calculation
#'   \item power spectrum density analysis including roughness against lengthscale calculation
#'   \item fractal dimension analysis including fractal dimensions calculation
#'   \item basic roughness parameters analysis such as mean, variance, Rrms, Ra
#' }
#'
#' @param AFMImageAnalyser a \code{\link{AFMImageAnalyser}} to manage and store image analysis
#' @return an \code{\link{AFMImageAnalyser}} containing all the analysis
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' AFMImage<-extractAFMImage(AFMImageOfAluminiumInterface, 0, 0, 32)
#' AFMImageAnalyser<-new("AFMImageAnalyser", AFMImage= AFMImage, fullfilename = AFMImage@@fullfilename)
#' AFMImageAnalyser<-analyse(AFMImageAnalyser)
#' print(AFMImageAnalyser@@fdAnalysis)
#' 
analyse<-function(AFMImageAnalyser) {
  AFMImage<-AFMImageAnalyser@AFMImage
  # Variogram analysis  
  #TODO 
  sampleFitPercentage<-3.43/100
  #sampleFitPercentage<-0.13/100
  variogramAnalysis<-AFMImageVariogramAnalysis(sampleFitPercentage)
  variogramAnalysis@omnidirectionalVariogram<- calculateOmnidirectionalVariogram(AFMImageVariogramAnalysis= variogramAnalysis, AFMImage=AFMImage)
  variogramAnalysis@directionalVariograms<- calculateDirectionalVariograms(AFMImageVariogramAnalysis= variogramAnalysis, AFMImage=AFMImage)
  
  # manage model evaluations
  AFMImageVariogram<-variogramAnalysis@omnidirectionalVariogram
  class(AFMImageVariogram)=c("gstatVariogram","data.frame")
  variogramAnalysis<-evaluateVariogramModels(variogramAnalysis, AFMImage)
  
  
  # PSD analysis
  psdAnalysis<-AFMImagePSDAnalysis()
  roughnessAgainstLengthscale(psdAnalysis)<-RoughnessByLengthScale(AFMImage, psdAnalysis)
  tryCatch({
    intersection <- getAutoIntersectionForRoughnessAgainstLengthscale(psdAnalysis, AFMImage, second_slope= FALSE)
    intersections<-c(intersection)
    intersection <- getAutoIntersectionForRoughnessAgainstLengthscale(psdAnalysis, AFMImage, second_slope= TRUE)
    intersections<-c(intersections,intersection)
    intersections(psdAnalysis)<-intersections
  }, error = function(e) {print(paste("Impossible to find PSD intersections automaticaly",e))})
  
  # fractal dimension analysis
  fdAnalysis<-AFMImageFractalDimensionsAnalysis()
  fractalDimensionMethods(fdAnalysis)<-getFractalDimensions(AFMImage, fdAnalysis)
  
  # basic roughness parameters
  AFMImageAnalyser@mean=mean(AFMImage@data$h)
  AFMImageAnalyser@variance=var(AFMImage@data$h)
  AFMImageAnalyser@TotalRrms=sqrt(var(AFMImage@data$h))
  AFMImageAnalyser@Ra=mean(abs(AFMImage@data$h))
  
  AFMImageAnalyser@variogramAnalysis<-variogramAnalysis
  AFMImageAnalyser@psdAnalysis<-psdAnalysis
  AFMImageAnalyser@fdAnalysis<-fdAnalysis
  
  return(AFMImageAnalyser)
}

#' Export all data from an analysis of an AFM Image as rdata file
#' 
#' A function to export to all the data from all analysis of an \code{\link{AFMImage}} and put them on disk as rdata file
#'
#' @param AFMImageAnalyser an \code{\link{AFMImageAnalyser}}
#' @param AFMImage an \code{\link{AFMImage}}
#' 
#' @name putAnalysisOnDisk
#' @rdname putAnalysisOnDisk-methods
#' @exportMethod putAnalysisOnDisk
#' @author M.Beauvais
#' 
setGeneric(name= "putAnalysisOnDisk", 
           def= function(AFMImageAnalyser, AFMImage) {
             return(standardGeneric("putAnalysisOnDisk"))
           })

#' @rdname putAnalysisOnDisk-methods
#' @aliases putAnalysisOnDisk,AFMImageAnalyser-method
setMethod(f="putAnalysisOnDisk", "AFMImageAnalyser",
          definition= function(AFMImageAnalyser, AFMImage) {
            filename<-basename(AFMImage@fullfilename)
            exportDirectory<-paste(dirname(AFMImage@fullfilename), "outputs", sep="/")
            
            saveAFMImageAnalyser(AFMImageAnalyser, AFMImage, exportDirectory)
            saveOnDisk(AFMImage, exportDirectory) # save AFMImage as rdata file
            
          })



setGeneric(name= "saveAFMImageAnalyser", 
           def= function(AFMImageAnalyser, AFMImage, exportDirectory) {
             return(standardGeneric("saveAFMImageAnalyser"))
           })

setMethod(f="saveAFMImageAnalyser", "AFMImageAnalyser",
          definition= function(AFMImageAnalyser, AFMImage, exportDirectory) {
            filename<-basename(AFMImage@fullfilename)
            
            exportCsvFilename<-paste(filename,"AFMImageAnalyser.rda", sep="-")
            exportCsvFullFilename<-paste(exportDirectory, exportCsvFilename, sep="/")
            print(paste("saving", basename(exportCsvFullFilename)))
            tryCatch({
              newVariableName=paste(filename,"AFMImageAnalyser.rda", sep="-")
              assign(newVariableName, AFMImageAnalyser)    
              save(list=c(newVariableName), file=exportCsvFullFilename)
            }, error = function(e) {print("error",e)})
            
          })


#' Calculate the total Root Mean Square Roughness (Rrms total)
#' 
#' \code{totalRMSRoughness} returns the total RMS roughness calculated from the variance of heights

#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @return a numeric as the square root of the variance of heights
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' totalRMSRoughness<-totalRMSRoughness(AFMImageOfAluminiumInterface)
#' print(totalRMSRoughness)
#' 
totalRMSRoughness<-function(AFMImage) {
  sqrt(var(AFMImage@data$h))
}

#' Get Roughness parameters
#'
#' Get basic roughness parameters as amplitude parameters:
#' Total root mean square Roughness or Total Rrms or totalRMSRoughness_TotalRrms\cr
#' Mean roughness or Ra or MeanRoughness_Ra
#'
#' \code{getRoughnessParameters} returns a data.table of roughness parameters
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @return a data.table of roughness parameters: 
#' \itemize{
#'   \item totalRMSRoughness_TotalRrms the total RMS Roughness as the square root of the variance of heights
#'   \item MeanRoughness_Ra the average roughness as the mean of absolute value of heights
#' }
#' @author M.Beauvais
#' @name getRoughnessParameters
#' @rdname getRoughnessParameters-methods
#' @exportMethod getRoughnessParameters
#' @examples
#' library(AFM)
#' 
#' data(AFMImageOfAluminiumInterface)
#' roughnessParameters<-getRoughnessParameters(AFMImageOfAluminiumInterface)
#' print(roughnessParameters)
#' 
setGeneric(name= "getRoughnessParameters", 
           def= function(AFMImage) {
             return(standardGeneric("getRoughnessParameters"))
           })

#' @rdname getRoughnessParameters-methods
#' @aliases getRoughnessParameters,AFMImage-method
setMethod(f="getRoughnessParameters", "AFMImage",
          definition= function(AFMImage) {
            # amplitude parameters
            totalRMSRoughness_TotalRrms = sqrt(var(AFMImage@data$h))
            MeanRoughness_Ra = mean(abs(AFMImage@data$h))
            #MeanRoughnessDepth_RzDIN = 
            #MaxProfileValleyDepth_Rmax=  
            
            # spacing parameters
            
            # hybrid parameters
            
            return(data.table(totalRMSRoughness_TotalRrms, MeanRoughness_Ra))
          })


#' Check the isotropy of a sample
#' 
#' \code{checkIsotropy} is used to check the isotropy of an \code{\link{AFMImage}}. 
#' A directional variogram is calculated for various directions. 
#' If the variogram is very similar for all the directions then the sample is isotropic.
#'
#' @param AFMImage an \code{\link{AFMImage}} to be analysed
#' @param AFMImageAnalyser an \code{\link{AFMImageAnalyser}} to perform the analysis
#' @return an \code{\link{AFMImageAnalyser}} containing the directional variograms
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' library(ggplot2)
#' 
#' data(AFMImageOfAluminiumInterface)
#' AFMImage<-extractAFMImage(AFMImageOfAluminiumInterface, 0, 0, 32)
#' AFMImageAnalyser<-new("AFMImageAnalyser", AFMImage= AFMImage, fullfilename = AFMImage@@fullfilename)
#' AFMImageAnalyser<-checkIsotropy(AFMImage,AFMImageAnalyser)
#' varios<-AFMImageAnalyser@@variogramAnalysis@@directionalVariograms
#' p2 <- ggplot(varios, aes(x=dist, y=gamma,  
#'                          color= as.factor(dir.hor), shape=as.factor(dir.hor)))
#' p2 <- p2 + expand_limits(y = 0)
#' p2 <- p2 + geom_point()
#' p2 <- p2 + geom_line()
#' p2 <- p2 + ylab("semivariance (nm^2)")
#' p2 <- p2 + xlab("distance (nm)")
#' p2 <- p2 + ggtitle("Directional")
#' p2
#' 
checkIsotropy<-function(AFMImage, AFMImageAnalyser) {
  print("checking isotropy...")
  
  # Variogram analysis  
  sampleFitPercentage<-3.43/100
  variogramAnalysis<-AFMImageVariogramAnalysis(sampleFitPercentage)
  if(!is.null(AFMImageAnalyser@updateProgress)) variogramAnalysis@updateProgress<-AFMImageAnalyser@updateProgress
  variogramAnalysis@directionalVariograms<- calculateDirectionalVariograms(AFMImage=AFMImage, AFMImageVariogramAnalysis=variogramAnalysis)
  
  AFMImageAnalyser@variogramAnalysis<-variogramAnalysis
  print("done.")
  return(AFMImageAnalyser)
}

#' Check visualy of the normality of the sample
#'
#' \code{checkNormality} performs a visual check to know if the distribution of heights of an \code{\link{AFMImage}} follows a normal distribution. The function displays Quantile/Quantile and distribution plots.
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param ... pngfullfilename (optional): directory and filename to save the visual check to png or pdffullfilename(optional): directory and filename to save the visual check to pdf
#' 
#' @references  Olea2006, Ricardo A. Olea "A six-step practical approach to semivariogram modeling", 2006, "Stochastic Environmental Research and Risk Assessment, Volume 20, Issue 5 , pp 307-318"
#'
#' @author M.Beauvais
#' @export
#' @examples
#' \dontrun{
#' library(AFM)
#' 
#' # display Quantile/Quantile and distribution plots.
#'   data(AFMImageOfNormallyDistributedHeights)
#'   checkNormality(AFMImage= AFMImageOfNormallyDistributedHeights)
#' 
#' # display and save on disk Quantile/Quantile and distribution plots.
#'   data(AFMImageOfNormallyDistributedHeights)
#'   checkNormality(AFMImage= AFMImageOfNormallyDistributedHeights, 
#'                  pngfullfilename=paste(tempdir(), "checkNormality.png", sep="/"))
#' }
checkNormality<- function(..., AFMImage) {
  force(AFMImage)
  
  fullfilename<-AFMImage@fullfilename
  
  args<-names(list(...))
  toPng<-c(match('pngfullfilename',args)!=-1)
  toPdf<-c(match('pdffullfilename',args)!=-1)
  
  qq <- checkNormalityQQ(AFMImage)
  m <- checkNormalityDensity(AFMImage)
  
  toFile<-FALSE
  if (!is.na(toPng)) {
    reportName<-paste(fullfilename, "-normality-checks",".png",sep="")
    print(paste("saving", basename(reportName)))
    png(reportName, width=1280, height=800, res=200)
    toFile<-TRUE
  }
  if (!is.na(toPdf)) {
    reportName<-paste(fullfilename, "-normality-checks",".pdf",sep="")
    print(paste("saving", basename(reportName)))
    pdf(reportName, width=11.69, height=8.27)
    toFile<-TRUE
  }
  
  if (toFile) {
    title<- paste("Normality tests for ",basename(fullfilename))
  }else{
    title<-"Normality tests"
  }
  
  grid.newpage() # Open a new page on grid device
  pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.5, 5, 0.5), "null"))))
  print(qq, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(m, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  other<-paste("mean=", mean(AFMImage@data$h), "- skewness=", moments::skewness(AFMImage@data$h))
  grid.text(other, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  
  if (toFile) {
    dev.off()
    #print("exist")
  }
}

checkNormalityQQ<- function(AFMImage) {
  h<-NULL
  ggplot(data=AFMImage@data, mapping=aes(sample=h)) + 
    stat_qq()  +  
    geom_abline(intercept = mean(AFMImage@data$h), slope = sd(AFMImage@data$h)) 
}

checkNormalityDensity<- function(AFMImage) {
  h<-..density..<-NULL
  ggplot(AFMImage@data, aes(x=h)) +
    geom_histogram( aes(y=..density..), colour="black", fill="white") +
    stat_function(fun = dnorm, args = list(mean = mean(AFMImage@data$h), sd =  sd(AFMImage@data$h)))
}

getLibrariesVersions<-function() {
  v<-data.table(installed.packages())
  Package<-NULL
  AFMPackage<-v[Package=="AFM",]
  gstatPackage<-v[Package=="gstat",]
  fractaldimPackage<-v[Package=="fractaldim",]
  fftwtoolsPackage<-v[Package=="fftwtools",]
  
  versions=data.table(lib=c("AFM",
                            "gstat",
                            "fractaldim", 
                            "fftwtools"),
                      version=c(AFMPackage$Version,
                                gstatPackage$Version,
                                fractaldimPackage$Version,
                                fftwtoolsPackage$Version)
                      )
  return(versions)
}
