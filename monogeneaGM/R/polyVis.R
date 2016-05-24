#' Wireframe plot of monogenean anchor
#'
#' This function plots the ventral and dorsal anchors as polygons in their natural positions on the specimen slide.  
#' @param index the index of a (tps) text file or matrix containing the (eleven) landmark coordinates in the following order: right ventral, left ventral, right dorsal, left dorsal
#' @param spacing a numeric constant specifying the spacing on the x and y coordinates relative to slide center
#' @param havelist choose TRUE if the landmark coordinate data are contained in a list
#' @param listdata a list containing objects that are matrices of 44 rows and 2 columns containing raw landmark coordinate data
#' @param tit title for the plot
#' @details This plot is useful for detecting slides with inconsistencies in magnification.  
#' It is also useful for detecting poor quality samples, as indicated by large shape variation between left and right forms. 
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' data(ligophorus_tpsdata)
#' #Simultaneously visualise the first four specimen slides
#' par(mfrow=c(2,2))
#' par(mar=c(3,3,2,1.5))
#' mapply(function(k) polyVis(k, spacing=300, havelist=TRUE, 
#' listdata=ligophorus_tpsdata$grandis), k=1:4)
#'

polyVis <- function(index, spacing=250, havelist=FALSE, listdata=NULL, tit=""){

if(havelist==TRUE){

x <- listdata[[index]]
centroid <- apply(x,2,mean)
plot((centroid[1]-spacing):(centroid[1]+spacing), main=tit,
(centroid[2]-spacing):(centroid[2]+spacing), pch="", xlab="x", ylab="y", xaxt="n", yaxt="n")

for(i in 1:4){
	polygon(x[(11*(i-1)+1):(11*i),1], x[(11*(i-1)+1):(11*i),2])
}

}

else {
x <- matrix(scan(dir()[index], skip=1, nlines=44),ncol=2,byrow=TRUE)
centroid <- apply(x,2,mean)

plot((centroid[1]-spacing):(centroid[1]+spacing), main=tit,
(centroid[2]-spacing):(centroid[2]+spacing), pch="", xlab="x", ylab="y", xaxt="n", yaxt="n")

for(i in 1:4){
	polygon(x[(11*(i-1)+1):(11*i),1], x[(11*(i-1)+1):(11*i),2])
}

}

}
