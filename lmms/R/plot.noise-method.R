# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the graphics and stats package.
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

#' Plot of \code{associations} objects
#' 
#' Plot of the filter ratios R_T and R_I as proposed by Straube et al 2014.
#' 
#' @import gplots
#' @import ggplot2
#' @param x an object of class \code{matrix} or \code{data.frame}.
#' @param colorBy \code{associations} the variable to be colored by.  Default  \code{'propMissing'}, options: \code{'propMissing'},\code{'fc'}.
#' @param fcCutoff an optional \code{numeric} value to remove ratios with low fold changes.
#' @param propMissingCutoff an optional \code{numeric} value to remove ratios with high number of missing values. 
#' @param \ldots ignored
#' @return plot showing filter ratios R_T and R_I as proposed by Straube \emph{et al.} 2014. Filter ratios can either be colored by proportion of missing values or maximum fold change. 
#' @references  Straube J., Gorse D., Huang B.E., Le Cao K.-A.(2014).  \emph{A linear mixed model spline framework for analyzing time course 'omics' data} Submitted
#' @seealso \code{\link{investNoise}}, \code{\link{filterNoise}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' G1 <- kidneySimTimeGroup$group=="G1"
#' noiseTest <-investNoise(data=kidneySimTimeGroup$data[G1,],time=kidneySimTimeGroup$time[G1],
#'             sampleID=kidneySimTimeGroup$sampleID[G1])
#' plot(noiseTest,colorBy="fc")
#' }
#' @method plot noise
#' @export
plot.noise <- function(x, colorBy="propMissing", fcCutoff=NA, propMissingCutoff=NA, ...){
indexFC = indexMissing = rep(TRUE,length(x@RT))

if(!is.na(fcCutoff)){
  rang <- range(x@foldChange,na.rm = T)
  if(!is.numeric(fcCutoff) | fcCutoff<rang[1] | fcCutoff>rang[2]){
    stop(paste("fcCutoff needs to be numeric and in range",rang[1],"-",rang[2]))
  }else{
    indexFC[x@foldChange<=fcCutoff] <- FALSE
  }
}

if(!is.na(propMissingCutoff)){
  rang <- range(x@propMissing,na.rm = T)
  if(!is.numeric(propMissingCutoff) | propMissingCutoff<rang[1] |propMissingCutoff>rang[2]){
    stop(paste("propMissingCutoff needs to be numeric and in range",rang[1],"-",rang[2]))
  }else{
    indexMissing[x@propMissing>propMissingCutoff] <- FALSE
  }
}
index <- indexFC & indexMissing
q <-qplot(x@RT[index],x@RI[index], col=color,xlab=expression(R[T]),ylab=expression(R[I])) + geom_point(alpha=0.5,size=3)
if(colorBy=="propMissing"){
  color <- x@propMissing[index]
  main = "Proportion of\nmissing\nvalues"
  q<-q +scale_colour_gradientn(limits=c(0,1),colours=c("#1b9e77","#d95f02","#7570b3"),name=main,breaks=seq(0,1, length.out = 5),labels=format(seq(0,1, length.out = 5)))  + theme_bw()
}else if(colorBy=="fc"){
  color <- x@foldChange[index]
  main = "Fold\nChange"
  minC <- min(color,na.rm=T)
  maxC <- max(color,na.rm=T)
  mc <- c()
  col <- c("black","green")
  if(minC<0){
    mc <- seq(minC,0,length.out = 3)
    col <-c("red","black","green")
  }
  rangeC <- c(mc,seq(0,maxC,length.out = 3))
  q <- q+scale_colour_gradientn(colours=col,name=main,breaks=rangeC,labels=format(signif(rangeC,2)))  + theme_bw()
}else{
  stop("Can not identify type of colorBy choose options: 'fc'or 'propMissing' ")
}
 q

}