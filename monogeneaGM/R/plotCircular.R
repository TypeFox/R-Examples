#' Circular plot
#'
#' This function creates a circular plot (Batschelet, 1981) showing how mean directional change of GPA coordinates in species
#' of interest is distributed relative the root ancestor. Magnitude of directional change is proportional to the length of
#' ray projecting from a data point.
#' @param x a list containing objects that are matrices of average GPA coordinates for a set of 
#' species. The names of this \code{x} contains should contain the species names
#' @param ancestor a matrix specifying the GPA coordinates of the root ancestor estimated
#' using the \code{fastAnc} function in the \code{phytools} package 
#' @param col.lab a character vector specifying colors for species with matching indices in \code{x}
#' @param clade a character vector specifying the species that form a clade of interest; currently supports only two clades
#' @param LM the landmark of interest
#' @param f scaling factors for the magnitude of directional change
#' @param sf shrinking factor for \code{shrink argument} in \code{plot.circular} from the \code{circular} package
#' @param ptscale scaling factor for size of data point on the circle perimeter
#' @param tit title for the circular plot
#' @details The arms in the circle are color-annotated according to clade of interest, with their direction and length indicating 
#' mean directional and magnitude change in the clades of interest, respectively. The circular plot is a useful complement to the wireframe-lollipop 
#' plot produced using \code{shapeEvo}, as it shows directional distribution details for all species at the individual landmark level. When
#' used in conjunction with the principal component loadings heat map produced using \code{pcloadhm}, it can greatly aid biological 
#' interpretation of the principal components. This function depends on the \code{phytools} (Revell, 2012) and \code{circular} (Agostinelli & Lund, 2013) packages. 
#' @seealso \code{\link{anglecheck}}, \code{\link{shapeEvo}}, \code{\link{pcloadhm}} 
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Agostinelli C, Lund U. (2013). R package 'circular': Circular Statistics (version 0.4-7). Available at: https://r-forge.r-project.org/projects/circular.
#'
#' Batschelet E. (1981). Circular Statistics in Biology. London: Academic Press.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217-223.
#' @examples
#' library(circular)
#'
#' data(va_mean)  
#' data(estimated_ancestral_va)
#' data(spcolmap)
#'
#' #the species in the defined clade infects the fish host Liza subviridis (dodger blue color)
#' plotCircular(va_mean, estimated_ancestral_va, col.lab=spcolmap$color, 
#' clade=spcolmap[spcolmap$host %in% "L.subviridis",]$species, LM=6, tit=6)
#'

plotCircular <- function(x, ancestor, col.lab, clade, LM, f=c(100,0.1), sf=1.5, ptscale=1, tit=NULL) {

#This is just to suppress annoying warning messages from functions in the circular package
#reminding users that the data is angular with radians unit

options(warn=-1)
dev <- lapply(x, function(k) k - ancestor)

angle_data <- lapply(dev, function(k) apply(k,1,anglecheck) )
angle_data <- do.call(rbind, angle_data)

cladeid <- which(names(x) %in% clade)

magnitude <- lapply(dev,function(k) apply(k, 1, function(h) sqrt ( sum(h^2) ) ) )
magnitude <- do.call(rbind, magnitude)

magnitude <- round( magnitude * f[1], 0) * f[2]

plot.circular(angle_data[,LM], shrink=sf, axes=FALSE, tol=0.01, main=tit)

cladeid1 <- which(names(dev) %in% clade)
cladeid2 <- setdiff(1:length(x), cladeid1)
for(i in cladeid1){
arrows.circular( angle_data[i,LM], col=col.lab[i], lwd=1, shrink=1+magnitude[i,LM], length=0 )
arrows.circular( angle_data[i,LM], col="white", lwd=2, length=0 )

points.circular( angle_data[i,LM], pch=16, col=col.lab[i], cex=1*ptscale)
points.circular( angle_data[i,LM], pch=1, cex=1*ptscale)

}

for(i in cladeid2){
arrows.circular( angle_data[i,LM], col=col.lab[i], lwd=1, shrink=1+magnitude[i,LM], length=0 )
arrows.circular( angle_data[i,LM], col="white", lwd=2, length=0 )

points.circular( angle_data[i,LM], pch=15, col=col.lab[i], cex=1*ptscale)
points.circular( angle_data[i,LM], pch=0, cex=1*ptscale)
}

text(0,0,"+")
text(0.8,0.0,0, cex=0.8)
text(0.0,0.7,expression(paste(italic(pi),"/2")),cex=0.8)
text(-0.8,0.0,expression(paste(italic(pi))),cex=0.8)
text(0.0,-0.7,expression(paste(3*italic(pi),"/2")),cex=0.8)

all.med <- mean(magnitude[cladeid1,LM])+mean(magnitude[cladeid2,LM])

arrows.circular(mean.circular(angle_data[cladeid1,LM]), length=0.05, lwd=1, shrink=1*mean(magnitude[cladeid1,LM])/all.med, col="dodgerblue")

arrows.circular(mean.circular(angle_data[cladeid2,LM]), length=0.05, lwd=1, shrink=1*mean(magnitude[cladeid2,LM])/all.med, col="violetred")

}

