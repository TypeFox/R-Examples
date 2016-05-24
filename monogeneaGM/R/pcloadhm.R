#' Heat map representation of principal component loadings of shape variables
#'
#' This function produces a graphical representation of selected principal component loadings of shape variables in the form
#' of a heat map.
#' @param x a matrix containing the loadings of shape variables (row) of each principal component (column) 
#' @param pc a constant specifying the principal component of interest
#' @param sgn the sign of the loadings; this value should follow the one used in \code{\link{pca2d}}
#' @param nrow a constant indicating the number of landmarks defined in the anchors; defaults to 11
#' @param color.code a character vector of hex color codes that define the color palette; if left undefined, defaults to
#' the red-white-blue palette
#' @param ylab y-axis title for the plot
#' @param xlab x-axis title for the plot
#' @param yaxis if TRUE, the y-axis values are labelled
#' @param tit title for the plot
#' @details The sign and magnitude of loadings of shape variables for a particular principal component 
#' is important for the latter's biological interpretation. A heat map representation is an alternative to 
#' the usual manner of presenting them in tabular form, and may be more effective for presentation purpose.
#' Reference to the circular plots (\code{plotCircular}) for each landmark and the PCA plots 
#' can be very useful in determining the biological interpretation of a particular principal component.
#' @seealso \code{\link{plotCircular}}, \code{\link{pca2d}}, \code{\link{colorBar}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' @examples
#' library(phytools)
#'
#' data(ligophorus_shape)
#' data(ligotree)
#' data(spcolmap)
#'
#' shapev <- pca2d(ligophorus_shape[,1:22], sgn=1, labcol=spcolmap$color, 
#' phylo=TRUE, phy=ligotree, genus="L. ",
#' bound.y = c(-0.1, 0.1), bound.x1 =c(-0.15,0.2), bound.x2=c(-0.15,0.2))
#'
#' fff <- c(0,1,1,2,2,3,3,0,4,4)
#' nf <- layout(matrix(c(rep(0,length(fff)),rep(fff,5),rep(0,length(fff))),
#' 7,length(fff),byrow=TRUE))
#' layout.show(nf)
#' par(mar=c(5,4,4,1))
#' #the loadings for the first three PC of shape variables of the ventral anchors
#' pcloadhm(shapev$variable,sgn=1,pc=1,yaxis=TRUE,tit="VPC1")
#' pcloadhm(shapev$variable,sgn=1,pc=2,ylab="", tit="VPC2")
#' pcloadhm(shapev$variable,sgn=1,pc=3,ylab="", tit="VPC3")
#'
#' #add a colorbar for completeness
#' par(mar=c(5,2,4,3))
#' colorBar(min=-1, max=1)
#'

pcloadhm <- function(x,pc=1,sgn=1,nrow=11,color.code=NULL,ylab="Landmark",
xlab="Coordinate", yaxis=FALSE, tit=NULL){

if(is.null(color.code)) {
rwb <- c("#99000D","#FB6A4A","white","#6BAED6","#084594")
color.code<- colorRampPalette(rwb, space="Lab")(101)
}

mat <- matrix(x[nrow(x):1,pc], nrow=nrow)
mat <- mat[,2:1] * sgn

image(1:ncol(mat), 1:nrow(mat), t(mat), col=color.code, xaxt="n", yaxt="n", 
ylab=ylab, xlab=xlab, main=tit)
axis(1, 1:2, c("x","y"), tick=FALSE)

if(yaxis==TRUE){
axis(2, nrow(mat):1, 1:nrow(mat), tick=FALSE, las=1)
}
lines(c(0.5,0.5),c(0.5,11.5))
lines(c(0.5,2.5),c(11.5,11.5))

}
