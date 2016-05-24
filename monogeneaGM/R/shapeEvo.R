#' Statistical test of deviation from directional uniformity and estimation of average magnitude of directional change
#'
#' This function performs the Rayleigh test (Batschelet, 1981) for detecting deviation from 
#' uniformity of directional change at each landmark. Additionally, it estimates the mean magnitude of directional change, 
#' and then summarizes the result graphically using a wireframe-lollipop plot (Klingenberg, 2013).
#' @param x a list of objects that are matrices containing average GPA coordinates of anchor landmarks. The   
#' species names should be the names of this list 
#' @param ancestor a matrix specifying the GPA coordinates of the root ancestor, estimated
#' using \code{fastAnc} function in the \code{phytools} package 
#' @param col.lab color for arrows in the wireframe-lollipop plot
#' @param coltones color tones for p-values; defaults to red-black-green spectrum
#' @param clade a character vector specifying the species that form a clade of interest
#' @param exfac an expansion factor for the magnitude of direction change
#' @param tit title for the wireframe-lollipop plot
#' @return A wireframe-lollipop plot and a list containing:
#' \item{magnitude}{a matrix of the mean magnitude of directional change (column) of each landmark for each species (row)}
#' \item{pvalue}{p-values for each landmark from the Rayleigh test}
#' @seealso \code{\link{plotCircular}}, \code{\link{anglecheck}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Batschelet E. (1981). Circular Statistics in Biology. London: Academic Press.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Klingenberg CP. (2013). Visualizations in geometric morphometrics: 
#' how to read and how to make graphs showing shape changes. Hystrix 24:15-24.
#' @examples
#' library(gplots)
#' library(circular)
#' 
#' data(va_mean)
#' data(estimated_ancestral_va)
#' data(spcolmap)
#'
#' nf <- layout(matrix(c(1,1,1,2),1, 4,byrow=TRUE))
#' layout.show(nf)
#' 
#' cladeII <- spcolmap$species[spcolmap$host %in% "L.subviridis"]
#' shapeEvo(va_mean, estimated_ancestral_va, col.lab="dodgerblue",
#' clade=cladeII, exfac=2, tit="Ventral anchors")
#' #Some journals want the title to be left-adjusted, so set tit="" and then:
#' #title("a)", adj=0)
#'
#' #Add a nice color bar
#' par(mar=c(5,6,4,2))
#' colorBar(redgreen(101),min=0, max=1, tit="p-value")
#'

shapeEvo <- function(x, ancestor, col.lab="black", coltones=redgreen(101), clade, exfac=1, tit=NULL) {
#This is just to suppress annoying warning messages from functions in the circular package
#reminding users that the data is angular with radians unit
options(warn=-1)
dev <- lapply(x, function(k) k - ancestor)

angle_data <- lapply(dev, function(k) apply(k,1,anglecheck) )
angle_data <- do.call(rbind, angle_data)

cladeid <- which(names(x) %in% clade)

pval <- round(apply(angle_data, 2, function(k) rayleigh.test(k[clade])$p.value ),3)
names(pval) <- mapply(function(k) paste(c("L",k), collapse=""), k=1:length(pval))

color.pval <- coltones[round(pval * 100, 0)]

magnitude <- lapply(dev,function(k) apply(k, 1, function(h) sqrt ( sum(h^2) ) ) )
magnitude <- do.call(rbind, magnitude)
colnames(magnitude) <- names(pval)


plot(ancestor, asp=1, main=tit,cex.axis=0.9,xlab="x",ylab="y", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) 
abline(h=0, v=0, lty=2)
points(ancestor, col=color.pval, pch=16, cex=1.3)
#points(ancestor,  pch=1, cex=1.3)
polygon(ancestor)

arrows.circular(apply(angle_data, 2, function(k) mean.circular(k[cladeid])), 
y=exfac*apply(magnitude,2,function(k) mean(k[cladeid])), 
x0=ancestor[,1], y0=ancestor[,2], col=col.lab, lwd=1, length=0.05) 
output <- list(magnitude,pval)

names(output) <- c("magnitude","pvalue")
return(output)

}

