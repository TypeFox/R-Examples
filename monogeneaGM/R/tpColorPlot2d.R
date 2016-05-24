#' Color-annotated scatter plot with transparency control
#'
#' This function enhances a typical two-dimensional scatter plot by enabling transparency control of color annotation.
#' In addition, it provides option for phylogenetic tree superimposition.
#' @param x a two-column matrix with rownames (usually, the species names)
#' @param labcol a character vector specifying species colors
#' @param xlab title for the x-axis
#' @param ylab title for the y-axis
#' @param tit title for the plot
#' @param circlesize a numeric vector that controls the centroid symbol size; defaults to 1 if \code{centroid = FALSE}
#' @param tpfac a numeric vector specifying the transparency level (0 to 255) for individual data points and the label mean
#' @param xbound range of values on the x-axis
#' @param ybound range of values on the y-axis 
#' @param centroid if TRUE, plots the centroid for each species
#' @param phylo if TRUE, coordinates of ancestral nodes from a supplied phylogeny (\code{phy}) are estimated using \code{fastAnc}
#' from the \code{phytools} package, and edges between nodes are joined according to the topology specified in \code{phy}
#' @param phy an object of class \code{phylo} from the \code{ape} package
#' @param pointscale a constant for controlling the size of the plotted ancestral nodes
#' @details Transparency control of color-annotated data points reduces visual saturation caused by the use of solid colours, 
#' thus allowing species centroids to be accentuated in the plot. In addition, if a user-supplied phylogeny is given,
#' it is superimposed onto the plot. This function depends on the \code{phytools} (Revell, 2012) and \code{ape} (Paradis et al., 2004)
#' packages.
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Paradis E, Claude J & Strimmer K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
#'
#' Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217-223.
#' @examples
#' data(pwed_pd)
#' data(spcolmap)
#' 
#' pwed_pd_list <- matrix2list.2(pwed_pd)
#' lm1 <- c("V1_3","V1_5")
#' 
#' #scatter plot of LM1-LM3 length against LM1-LM5 length for the ventral anchors
#' tpColorPlot2d(pwed_pd[,colnames(pwed_pd) %in% lm1], labcol=spcolmap$color, 
#' xlab=expression(paste("V1_3 ", "(", italic(mu),"m", ")")),
#' ylab=expression(paste("V1_5 ", "(", italic(mu),"m", ")")), centroid=TRUE)
#'

tpColorPlot2d <- function (x, labcol = "", xlab = "", ylab = "", tit = "", 
    circlesize = rep(1.5,nrow(x)), tpfac = c(60, 255), xbound=NULL, ybound=NULL, 
    centroid=FALSE,  phylo=FALSE, phy, pointscale=1) {

    label <- levels(as.factor(rownames(x)))     
   
    if(is.null(xbound) == TRUE & is.null(ybound) == TRUE){

    xbound <- c(floor(min(x[, 1])), ceiling(max(x[, 1])))
    ybound <- c(floor(min(x[, 2])), ceiling(max(x[, 2])))

    }

    plot(x[, 1], x[, 2], cex = circlesize[1], col = "white", 
        main = tit, xlab = xlab, ylab = ylab, xlim = xbound, 
        ylim = ybound)

	if(phylo==TRUE){

	mean.mat <- matrix(0,length(label),ncol(x))

	for(i in 1:ncol(x)){
	temp <- stack(x[,i])
	mean.mat[,i] <- tapply(temp[,1], temp[,2], mean)
	}

	rownames(mean.mat) <- label
	#ancestral node estimation
	#only the first three coordinates considered, unless more are requested by user

	anc.coord <- matrix(0,phy$Nnode,2)

	for(i in 1:2){
	anc.coord[,i] <- fastAnc(phy,mean.mat[,i])
	}

	#do lines first
	all.nodes <- mean.mat
	all.nodes <- all.nodes[phy$tip.label,]

	all.nodes <- rbind(all.nodes, anc.coord)

	for (i in 1:nrow(phy$edge)) {
        lines(all.nodes[phy$edge[i,], 1], all.nodes[phy$edge[i,], 2])
      }

	for(i in 1:nrow(anc.coord)){
	points(anc.coord[i,1], anc.coord[i,2], pch=16, cex=0.5 * pointscale, col="white")
	points(anc.coord[i,1], anc.coord[i,2], cex=0.5 * pointscale)
	}

	}

    for (i in 1:nrow(x)) {
        sp <- which(rownames(x) == label[i])
        red.green.blue <- as.numeric(col2rgb(labcol[i]))
        points(x[sp, 1], x[sp, 2], col = rgb(red = red.green.blue[1], 
            green = red.green.blue[2], blue = red.green.blue[3], 
            alpha = tpfac[1], maxColorValue = 255), pch = 16)

	   if(centroid==TRUE){
        points(mean(x[sp, 1]), mean(x[sp, 2]), col = rgb(red = red.green.blue[1], 
            green = red.green.blue[2], blue = red.green.blue[3], 
            alpha = tpfac[2], maxColorValue = 255), pch = 16, cex = circlesize[i])
        points(mean(x[sp, 1]), mean(x[sp, 2]), pch = 1, cex = circlesize[i])
	   }
    }
   

}

