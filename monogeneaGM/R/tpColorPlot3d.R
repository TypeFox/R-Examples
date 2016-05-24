#' Color-annotated three-dimensional scatter plot with transparency control
#'
#' This function enhances a three-dimensional scatter plot by enabling transparency control of 
#' color annotation. In addition, it provides option for phylogenetic tree superimposition.
#' @param x a matrix with rows representing samples and columns representing three variables of
#' interest, typically principal components 
#' @param r radius of plotting sphere 
#' @param phylo if TRUE, coordinates of ancestral nodes from a supplied phylogeny (\code{phy}) are estimated using \code{fastAnc}
#' from the \code{phytools} package, and edges between nodes are joined according to the topology specified in \code{phy}
#' @param phy an object of class \code{phylo} from the \code{ape} package
#' @param labcol a character vector specifying species colors
#' @param xyzlabel a vector of characters specifying the titles for the xyz-axes
#' @param alpha.set a constant for controlling degree of transparency (0 for complete transparency; 1 for solid color) of the data points
#' @param mean.show if TRUE, the centroids of each species is plotted in solid color
#' @param asp a vector specifying the aspect ratio of the xyz axes; the default gives a cube
#' @details Transparency control of color-annotated data points reduces visual saturation caused by the use of solid colors for all 
#' data points. Guide for choosing optimal value of \code{r}: for data range between -0.1 and 0.1, a value of 0.005 should be adequate. 
#' If a phylogenetic tree is supplied, it may be superimposed onto the three-dimensional space to allow visualization of 
#' evolutionary trajectories.
#' @seealso \code{\link{tpColorPlot2d}}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Paradis E, Claude J & Strimmer K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
#'
#' Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217-223.
#' @examples
#' library(phytools)
#' library(rgl)
#'
#' data(ligophorus_shape)
#' data(ligotree)
#' data(spcolmap)
#'
#' #Perform PCA of shape data for dorsal anchors and make 2D plots
#' pcashape <- pca2d(ligophorus_shape[,23:44], labcol=spcolmap$color, 
#' phylo=TRUE, phy=ligotree, genus="L. ", bound.y=c(-0.08, 0.1), 
#' bound.x1=c(-0.2,0.2), bound.x2 = c(-0.2,0.2))
#'
#' #Check for proportion of variation explained by each PC
#' summary(pcashape$pca)
#'
#' #A closer look with 3D plot
#' tpColorPlot3d(pcashape$scores[,3:1], r=0.005, phylo=TRUE, phy=ligotree, labcol=spcolmap$color,
#' xyzlabel=c("PC3 (8%)","PC2 (10%)","PC1 (61%)"), mean.show=TRUE)
#'

tpColorPlot3d <- function(x,r=0.005,phylo=FALSE,phy,labcol,xyzlabel=NULL,alpha.set=0.2, 
mean.show=FALSE,asp=c(1,1,1)) {
sp <- levels(as.factor(rownames(x)))

mean.mat <- matrix(0,length(sp),ncol(x))

for(i in 1:ncol(x)){
	temp <- stack(x[,i])
	mean.mat[,i] <- tapply(temp[,1], temp[,2], mean)
}
rownames(mean.mat) <- sp

for(i in 1:length(sp)){
spheres3d(x[rownames(x) %in% sp[i],1], x[rownames(x) %in% sp[i],2], 
x[rownames(x) %in% sp[i],3], col=labcol[i], type="s", radius=r, add=TRUE, alpha=alpha.set)
}
aspect3d(asp[1],asp[2],asp[3])

if(phylo == TRUE){
anc.coord <- matrix(0,phy$Nnode,3)

for(i in 1:3){
	anc.coord[,i] <- fastAnc(phy,mean.mat[,i])
	}

#do lines first
all.nodes <- mean.mat
all.nodes <- all.nodes[phy$tip.label,]

all.nodes <- rbind(all.nodes, anc.coord)

for (i in 1:nrow(phy$edge)) {
        lines3d(all.nodes[phy$edge[i, ], 1], all.nodes[phy$edge[i, ], 2], all.nodes[phy$edge[i, ], 3], lwd=1)
    }

}

if(mean.show==TRUE){
for(i in 1:length(sp)){
spheres3d(mean.mat[rownames(mean.mat) %in% sp[i],1], mean.mat[rownames(mean.mat) %in% sp[i],2], 
mean.mat[rownames(mean.mat) %in% sp[i],3], col=labcol[i], type="s", radius=r*1.5, add=TRUE)
}
}

rgl.bbox(color="grey50", emission="grey50", xlen=0, ylen=0, zlen=0)
rgl.material(color="black")

axes3d(edges=c("x--","y+-","z--"), nticks=4, cex=0.75)

mtext3d(xyzlabel[1], edge="x--", line=2)
mtext3d(xyzlabel[2], edge="y+-", line=3)
mtext3d(xyzlabel[3], edge="z--", line=3)

}