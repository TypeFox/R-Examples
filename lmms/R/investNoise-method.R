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

#' Quality control for time course profiles
#' 
#' Function to calculate filter ratios of trajectories.
#' 
#' @importFrom parallel parLapply detectCores makeCluster clusterExport stopCluster
#' @import methods
#' @importFrom stats sd
#' @usage investNoise(data, time, sampleID, log, numCores)
#' @param data \code{data.frame} or \code{matrix} containing the samples as rows and features as columns
#' @param time \code{numeric} vector containing the sample time point information.
#' @param sampleID \code{character}, \code{numeric} or \code{factor} vector containing information about the unique identity of each sample
#' @param log \code{logical} indicating log transformation of the data. Default value is TRUE
#' @param numCores alternative \code{numeric} value indicating the number of CPU cores to be used for parallelization. Default value is automatically estimated.
#' @details
#' investNoise calculates filter ratios R_T and R_I based on the time, individual and overall standard deviation as proposed by Straube \emph{et al.}  2015. 
#' @return investNoise returns an object of class \code{noise} containing the following components:
#' \item{name}{\code{character} the colnames or the index.}
#' \item{RT}{\code{numeric} the time to molecule sd ratio of each trajectory.} 
#' \item{RI}{\code{numeric} the individual to molecule sd ratio of each trajectory.}
#' \item{propMissing}{\code{numeric} Proportion of missing values for each trajectory. }
#' \item{foldChange}{\code{numeric} the maximum absolute fold change (either for log transformed data max(time)-min(time) or not log transformed data max(time)/min(time)) observed between the mean of any two time points. }
#' @references  Straube J., Gorse D., Huang B.E., Le Cao K.-A. (2015).  \emph{A linear mixed model spline framework for analyzing time course 'omics' data} PLOSONE, 10(8), e0134540.
#' @seealso \code{\link{summary.noise}}, \code{\link{plot.noise}}, \code{\link{filterNoise}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' G1 <- kidneySimTimeGroup$group=="G1"
#' noiseTest <-investNoise(data=kidneySimTimeGroup$data[G1,],time=kidneySimTimeGroup$time[G1],
#'             sampleID=kidneySimTimeGroup$sampleID[G1])
#' summary(noiseTest)
#' plot(noiseTest,colorBy="propMissing")}
# @docType methods
# @rdname investNoise-methods
# @export
# setGeneric('investNoise',function(data,time,sampleID,log,numCores){standardGeneric('investNoise')})
# setClassUnion("missingOrnumeric", c("missing", "numeric"))
# setClassUnion("missingOrlogical", c("missing", "logical"))
# setClassUnion("factorOrcharacterOrnumeric", c( "factor","character","numeric"))
# setClassUnion("matrixOrframe",c('matrix','data.frame'))
# # @rdname investNoise-methods
# # @aliases investNoise,matrixOrframe,numeric,factorOrcharacterOrnumeric,missingOrnumeric-method
# # @exportMethod investNoise
# 
# setMethod('investNoise',c(data="matrixOrframe",time="numeric",sampleID="factorOrcharacterOrnumeric",log='missingOrlogical',numCores="missingOrnumeric"), function(data,time,sampleID,log,numCores){
#   invest.Noise(data=data,time=time,sampleID=sampleID,log=log,numCores=numCores)
# })
#' @docType methods
#' @rdname investNoise-methods
#' @export
investNoise <- function(data,time,sampleID,log,numCores){
  nr <- nrow(data)
  nc <- ncol(data)
  
  if((nr != length(time)) & (nr != length(sampleID))){
    stop('Number of rows must be equal to the length of time and sampleID')
  }
  
  if(is.null(colnames(data))){
    warning("Colnames are missing. Varibale name is the column index.")
    name <- paste("X",1:nc,sep="")
  }else{
    name <- colnames(data)
  }
  if(missing(log))
    log <- T
    
  if(missing(numCores)){
    num.Cores <- detectCores()
  }else{
    num.Cores <- detectCores()
    if(num.Cores<numCores){
      warning(paste('The number of cores is bigger than the number of detected cores. Using the number of detected cores',num.Cores,'instead.'))
    }else{
      num.Cores <- numCores
    }
  }
    
  cl <- makeCluster(num.Cores,"SOCK")
  
  clusterExport(cl, list('data','time','sampleID'),envir=environment())
  
  new.data <- parLapply(cl,1:nc,fun = function(i){
    
  d <- data[,i]
  
  
  sd.mol <-sd(d,na.rm=T)
  mean.time <-tapply(d,time,function(x)mean(x,na.rm=T))
  if(log){
   fc <- abs(max(mean.time,na.rm=T)-min(mean.time,na.rm=T))
  }else{
    fc <- abs(max(mean.time,na.rm=T)/min(mean.time,na.rm=T))
  }
  
  sd.time <- mean(tapply(d,time,function(x)sd(x,na.rm=T)),na.rm=T)/sd.mol
  sd.ind <- 1-(mean(tapply(d,sampleID,function(x)sd(x,na.rm=T)),na.rm=T)/sd.mol)

  num.na <- sum(is.na(d))

  return(list(sd.time=sd.time,sd.ind=sd.ind,num.na=num.na,fc=fc))
  
  })
  
  stopCluster(cl)
  ratio.sd <- unlist(sapply(new.data,'[','sd.time'))
  ratioInd.sd <- unlist(sapply(new.data,'[','sd.ind'))
  tl <-length(unique(time))
  sl <- length(unique(sampleID))
  num.na <- unlist(sapply(new.data,'[','num.na'))/(tl*sl)
  fc <-  unlist(sapply(new.data,'[','fc'))

  l <- new('noise',name=name,RT =ratio.sd,RI=ratioInd.sd,propMissing=num.na, foldChange=fc)
  
  return(l)
}
