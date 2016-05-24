# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#' Filter non-informative trajectories
#' 
#' Function to remove non-informative trajectories
#' 
#' 
#' @usage filterNoise(data, noise, RTCutoff, RICutoff, propMissingCutoff, fcCutoff)
#' @param data \code{data.frame} or \code{matrix} containing the samples as rows and features as columns.
#' @param noise an object of class \code{noise} containing time and individual to molecule sd ratios number of missing values and maximum fold changes. 
#' @param RTCutoff \code{numeric} the R_T cutoff to remove non-informative trajectories.
#' @param RICutoff \code{numeric} the R_I to remove non-informative trajectories.
#' @param propMissingCutoff \code{numeric} maximum proportion of missing values in trajectories.
#' @param fcCutoff \code{numeric} the minimum fold change observed between the mean of any two time points. 
#' @details
#' filterNoise removes noisy or non-informative profiles based on selected theresholds R_I, R_T (Straube \emph{et al.} 2015), maximum foldchanges and/or missing values.
#' @return filterNoise returns an object of class \code{list} containing the following components:
#' \item{data}{\code{numeric} filtered data.}
#' \item{removedIndices}{\code{numeric} removed indices}
#' @references  Straube J., Gorse A.-D., Huang B.E., Le Cao K.-A. (2015).  \emph{A linear mixed model spline framework for analyzing time course 'omics' data} PLOSONE, 10(8), e0134540.
#' @seealso \code{\link{investNoise}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' G1 <- kidneySimTimeGroup$group=="G1"
#' noiseTest <-investNoise(data=kidneySimTimeGroup$data[G1,],time=kidneySimTimeGroup$time[G1],
#'             sampleID=kidneySimTimeGroup$sampleID[G1])
#' data <-filterNoise(data=kidneySimTimeGroup$data[G1,],noise=noiseTest,RTCutoff=0.9,
#'               RICutoff=0.3,propMissingCutoff=0.5)$data
#'              
#'              
#' #Alternatively model-based clustering can be used for filtering
#' library(mclust)
#' clusterFilter <- Mclust(cbind(noiseTest@@RT,noiseTest@@RI),G=2)
#' plot(clusterFilter,what = "classification")
#' meanRTCluster <-tapply(noiseTest@@RT,clusterFilter$classification,mean)
#' bestCluster <- names(meanRTCluster[which.min(meanRTCluster)])
#' filterdata <- kidneySimTimeGroup$data[G1,clusterFilter$classification==bestCluster]
#'               
#' }
#' @docType methods
#' @rdname filterNoise-methods
#' @export
setGeneric('filterNoise',function(data,noise,RTCutoff,RICutoff,propMissingCutoff,fcCutoff){standardGeneric('filterNoise')})

setClassUnion("matrixOrframe",c('matrix','data.frame'))
setClassUnion("missingOrnumeric", c("missing", "numeric"))
#' @rdname filterNoise-methods
#' @aliases filterNoise,matrixOrframe,noise,missingOrnumeric,missingOrnumeric,missingOrnumeric,missingOrnumeric-method
#' @exportMethod filterNoise

setMethod('filterNoise',c(data="matrixOrframe",noise="noise",RTCutoff="missingOrnumeric",RICutoff="missingOrnumeric",propMissingCutoff="missingOrnumeric",fcCutoff="missingOrnumeric"), function(data,noise,RTCutoff,RICutoff,propMissingCutoff,fcCutoff){
  filter.Noise(data=data,noise=noise,RTCutoff=RTCutoff,RICutoff=RICutoff,propMissingCutoff=propMissingCutoff,fcCutoff=fcCutoff)
})



filter.Noise <- function(data,noise,RTCutoff,RICutoff,propMissingCutoff,fcCutoff){
    rows <- ncol(data)
    justna <- colSums(is.na(data))!=rows
    
  if(all(missing(RTCutoff),missing(RICutoff),missing(propMissingCutoff),missing(fcCutoff))){
    RTCutoff=0.9
    RICutoff=0.3
  
    interIndex <- (noise@RT<=RTCutoff & noise@RI<=RICutoff)& justna
    
    l <- list(data=data[,interIndex],removedIndices=which(interIndex==F))
    return(l)
  }else{
    indexTime <- indexInd<-indexMissing <- indexFC <- rep(TRUE,dim(data)[2])
    if(!missing(propMissingCutoff)){
      indexMissing[noise@propMissing>=propMissingCutoff] <- FALSE
    }
    if(!missing(fcCutoff)){
      indexFC[noise@foldChange<=fcCutoff] <- FALSE
    }
    if(!missing(RICutoff)){
      indexInd[noise@RI>=RICutoff] <- FALSE
    }
    
    if(!missing(RTCutoff)){
      indexTime[noise@RT>=RTCutoff] <- FALSE
    }
  }
  
  interIndex <- (indexInd==indexTime) & indexFC & indexMissing & justna
  l <- list(data=data[,interIndex],removedIndices=which(interIndex==F))
  return(l)
  
}
