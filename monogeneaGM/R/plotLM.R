#' Scatter plot of Generalized Procrustes Analysis (GPA) coordinates of anchor landmarks
#'
#' This function generates a scatter plot of GPA coordinates of anchor landmarks, with the option of making a 
#' wireframe plot by joining mean GPA coordinates of adjacent landmarks. 
#' @param x an array containing GPA coordinate data of anchor landmarks 
#' @param tit title of the plot
#' @param pointscale a numeric constant for controlling the symbol size for observations
#' @param axispointscale a numeric constant for controlling the font size of labels of values on the xy axes
#' @param meansize a constant for controlling the symbol size of the mean coordinates
#' @param polygon.outline if TRUE, a wireframe plot connecting all adjacent landmarks is made
#' @param xbound a numeric vector specifying the range of x values for the plot
#' @param ybound a numeric vector specifying the range of y values for the plot
#' @param pch.suppress if TRUE, only the mean coordinates of the landmarks are plotted
#' @details The resulting scatter plot is an important graphical sanity check for potential problems after performing GPA.
#' @seealso \code{\link{procrustesFit}}, \code{\link{stdLM}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @examples
#' data(liewi_gpa)
#'
#' nice_title <- expression(paste("Dorsal anchor ", italic(L.liewi)))
#' plotLM(liewi_gpa, tit=nice_title, pointscale=0.8, axispointscale=0.8, 
#' meansize=1.2, polygon.outline=TRUE,c(-.6,.6),c(-.6,.6) )
#'

plotLM <- function(x, tit="", pointscale=1, axispointscale=1, 
meansize=1,polygon.outline=FALSE,xbound=NULL,ybound=NULL, pch.suppress=FALSE){

if(is.null(xbound) & is.null(ybound)){

if(pch.suppress==FALSE){
plot(x[,1,], x[,2,], pch=21, bg="gray", cex.axis=axispointscale,
xlab="x",ylab="y", cex=pointscale * 1, asp=1, main=tit)
}

else if(pch.suppress==TRUE){
plot(x[,1,], x[,2,], pch="", bg="gray", cex.axis=axispointscale,
xlab="x",ylab="y", cex=pointscale * 1, asp=1, main=tit)
}


}

else if(is.null(xbound) == FALSE | is.null(ybound) == FALSE){

if(pch.suppress==FALSE){
plot(x[,1,], x[,2,], pch=21, bg="gray", xlim=xbound, ylim=ybound,
xlab="x",ylab="y", cex=pointscale * 1, main=tit, cex.axis=axispointscale,)
}

else if(pch.suppress==TRUE){
plot(x[,1,], x[,2,], pch="", bg="gray", xlim=xbound, ylim=ybound,
xlab="x",ylab="y", cex=pointscale * 1, main=tit, cex.axis=axispointscale,)
}


}
abline(h=0, v=0, lty=2)

centroid <- apply(x, c(1,2), mean)
points(centroid, pch=21, bg="black", cex=meansize)

if(polygon.outline == TRUE){

polygon(centroid,lwd=1,border="black")

}

}