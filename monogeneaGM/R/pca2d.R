#' Principal component analyis 
#'
#' This function performs principal component analysis (PCA) and produces color-annotated scatter plots of a reference
#' principal component (typically the first one) against two other principal components (typically the second and third).
#' @param x a matrix with rows representing samples and columns representing morphometrical variables of
#' interest
#' @param a a vector of length 3 for the principal component ranks specified by the user; defaults to the first three principal components
#' @param sgn a numeric constant, either -1 or 1, that controls the sign of the principal component scores; defaults to 1
#' @param labcol a character vector giving the color-annotation of the species 
#' @param bound.x1 a numeric vector specifying the range of values on the x-axis for the first plot
#' @param bound.x2 a numeric vector specifying the range of values on the x-axis for the second plot
#' @param bound.y a numeric vector specifying the range of values on the y-axis for both plots
#' @param pointscale a constant for the size of species centroids; defaults to 1
#' @param phylo if TRUE, coordinates of ancestral nodes from a supplied phylogeny (\code{phy}) are estimated using \code{fastAnc} from the 
#' \code{phytools} package (Revell, 2012), and edges between nodes are joined according to the tree topology specified in \code{phy}
#' @param phy an object of class \code{phylo}
#' @param genus single character abbreviation for genus
#' @details To be specific, this function implements an R-mode, covariance-based PCA. When variables differ in their
#' units of measurement or show large magnitude differences, a correlation-based PCA is more reasonable. In this case,
#' the data in the input matrix must first be normalized by subtracting mean and dividing by standard deviation.
#' In R-mode PCA, the principal components can be interpreted contextually by checking the loadings of the variables of interest (graphically
#' using \code{pcloadhm}). However, this requires that the number of rows exceeds the number of columns (variables). If specimen sample size is small,
#' Q-mode PCA is possible. In this case, the input data matrix is transposed. Although species clusters can still be visualized, 
#' the principal components do not seem to be interpretable. If a phylogeny of the species is available, it can be superimposed onto the
#' principal component space to yield a phylomorphospace to provide a graphical complement to formal phylogenetic signal testing.
#' For the latter, see \code{physignal} in the \code{geomorph} package (Adams & Otarola-Castillo, 2013).
#' @seealso \code{\link{pcloadhm}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Adams DC, Otarola-Castillo E. (2013). geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393-399.
#'
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217-223.
#' @examples
#' library(phytools)
#'
#' data(ligotree)	
#' data(ligophorus_shape)
#' data(spcolmap)
#'
#' #PCA plot for the shape variables of the ventral anchors
#' pca2d(ligophorus_shape[,1:22], labcol=spcolmap$color, phylo=TRUE,
#' phy=ligotree, genus="L. ", bound.y=c(-0.1, 0.1), bound.x1=c(-0.2,0.2), 
#' bound.x2 = c(-0.2,0.2))
#'

pca2d <- function(x, a=1:3, sgn=1, labcol,
bound.x1=c(-0.15,0.15), bound.x2=c(-0.15,0.15), bound.y=c(-0.15,0.15),
pointscale=1, phylo=FALSE, phy, genus="") {

	a1 <- a[1] ; a2 <- a[2]; a3 <- a[3]
	splabel <- as.character(levels(as.factor(rownames(x))))
	nsp <- length(splabel)

	#Principal Component Analysis 
	pca <- princomp(x, cor=FALSE)
	lambda <- pca$sdev * sqrt(pca$n.obs)
	scores <- sgn * t ( t(pca$scores) / lambda )
	variable <- sgn * t ( t(pca$loadings) * lambda )

	nf <- layout(matrix(c(1,2,1,2,3,3),3, 2,byrow=TRUE))
	layout.show(nf)

	varexp <- round ( (pca$sdev)^2 / sum( (pca$sdev)^2 ) * 100, 0)

	par(mar=c(5,4,4,0) + 0.2)

	plot(scores[,a2], scores[,a1], xlab=paste("PC",a2,"(",varexp[a2],"%)"), 
	ylab=paste("PC",a1,"(",varexp[a1],"%)"), cex=1, pch="", xlim=bound.x1, ylim=bound.y)

	if(phylo==TRUE){

	nsp <- length(levels(as.factor(rownames(x))))

	mean.mat <- matrix(0,nsp,ncol(x))

	for(i in 1:ncol(x)){
		temp <- stack(scores[,i])
		mean.mat[,i] <- tapply(temp[,1], temp[,2], mean)
	}

	rownames(mean.mat) <- splabel

#ancestral node estimation
#only the first three coordinates considered, unless more are requested by user

	varname <- c(a1,a2)
	anc.coord <- matrix(0,phy$Nnode,2)

	for(i in c(1,2)){
		anc.coord[,i] <- fastAnc(phy,mean.mat[,varname[i]])
	}

#do lines first
	all.nodes <- mean.mat[,varname]
	all.nodes <- all.nodes[phy$tip.label,]

	all.nodes <- rbind(all.nodes, anc.coord)

	for(i in 1:nrow(phy$edge)) {
      	lines(all.nodes[phy$edge[i,], 2], all.nodes[phy$edge[i,], 1])
    	}


	for(i in 1:nrow(anc.coord)){
		points(anc.coord[i,2], anc.coord[i,1], pch=16, cex=0.8, col="white")
		points(anc.coord[i,2], anc.coord[i,1], cex=0.8)
	}

}

for(i in 1:nrow(scores)){
	sp <- which(rownames(scores)==splabel[i])
	red.green.blue <- as.numeric(col2rgb(labcol[i]))
	points(scores[sp,a2],scores[sp,a1],
	col=rgb(red=red.green.blue[1], green=red.green.blue[2], blue=red.green.blue[3], 
	alpha=60, maxColorValue=255), pch=16, cex=1)

	points(mean(scores[sp,a2]), mean(scores[sp,a1]), 
	col=rgb(red=red.green.blue[1], green=red.green.blue[2], blue=red.green.blue[3], 
	alpha=220, maxColorValue=255), pch=16, cex=1.5) 
	points(mean(scores[sp,a2]), mean(scores[sp,a1]), pch=1, cex=1.5)
}

par(mar=c(5,2,4,2) + 0.2)
plot(scores[,a3], scores[,a1], xlab=paste("PC",a3,"(",varexp[a3],"%)"), 
cex=1, yaxt="n", pch="", xlim=bound.x2, ylim=bound.y)


	if(phylo==TRUE){

	mean.mat <- matrix(0,nsp,ncol(x))

	for(i in 1:ncol(x)){
	temp <- stack(scores[,i])
	mean.mat[,i] <- tapply(temp[,1], temp[,2], mean)
	}

rownames(mean.mat) <- splabel
#ancestral node estimation
#only the first three coordinates considered, unless more are requested by user


anc.coord <- matrix(0,phy$Nnode,2)
varname <- c(a1,a3)
for(i in 1:2){
	anc.coord[,i] <- fastAnc(phy,mean.mat[,varname[i]])
	}

#do lines first
all.nodes <- mean.mat[,varname]
all.nodes <- all.nodes[phy$tip.label,]

all.nodes <- rbind(all.nodes, anc.coord)

for (i in 1:nrow(phy$edge)) {
        lines(all.nodes[phy$edge[i,], 2], all.nodes[phy$edge[i,], 1])
    }

for(i in 1:nrow(anc.coord)){
	points(anc.coord[i,2], anc.coord[i,1], pch=16, cex=0.8, col="white")
	points(anc.coord[i,2], anc.coord[i,1], cex=0.8)
}

}


for(i in 1:nrow(scores)){
	sp <- which(rownames(scores)==splabel[i])
	red.green.blue <- as.numeric(col2rgb(labcol[i]))
	points(scores[sp,a3],scores[sp,a1],
	col=rgb(red=red.green.blue[1], green=red.green.blue[2], blue=red.green.blue[3], 
	alpha=60, maxColorValue=255), pch=16, cex=1)

	points(mean(scores[sp,a3]), mean(scores[sp,a1]), 
	col=rgb(red=red.green.blue[1], green=red.green.blue[2], blue=red.green.blue[3], 
	alpha=255, maxColorValue=255), pch=16, cex=1.5) 
	points(mean(scores[sp,a3]), mean(scores[sp,a1]), pch=1, cex=1.5)
}


splabel_full <- mapply(function(k) paste(c(genus,splabel[k]), collapse=""), k=1:length(splabel)) 

f <- vector("expression",length(splabel_full))
for(s in 1:length(splabel_full)){
f[[s]] <- substitute(italic(nn), list(nn=splabel_full[s]))
}


par(mar=c(5,6,3,1))
plot(1:10,1:10,pch="", bty="n", xaxt="n", yaxt="n", main="Species", xlab="", ylab="")
legend(1,10,pch=rep(16,length(splabel)), pt.cex=2, col=labcol,f, ncol=4)
legend(1,10,pch=rep(1,length(splabel)), pt.cex=2,f, ncol=4)

output <- list(pca,lambda,scores,variable)
names(output) <- c("pca","lambda","scores","variable")

return(output) 
}

