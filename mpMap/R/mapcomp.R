#' Functions for comparison of two map orders
#' 
#' Takes in two maps with the aim of comparing the position of common markers. Creates a mapcomp object for plotting and summary. 
#' @rdname mapcomp-all
#' @aliases mapcomp plot.mapcomp summary.mapcomp
#' @usage mapcomp(object1, object2)
#' \method{plot}{mapcomp}(x, ...)
#' \method{summary}{mapcomp}(object, ...)
#' @export mapcomp
#' @param object1 Object inheriting class \code{mpcross} or class \code{map} 
#' @param object2 Object inheriting class \code{mpcross} or class \code{map}
#' @param object Object of class \code{mapcomp} for summarizing
#' @param x Object of class \code{mapcomp}
#' @param ... Additional arguments
#' @return An object of class \code{mapcomp} with components:
#' \item{commonmrk}{A matrix containing 5 columns with the names of all common markers for the two maps (mname), the chromosome mapped to in the first map (chr1), the position mapped to in the first map (pos1), the chromosome mapped to in the second map (chr2) and the position mapped to in the second map (pos2)}
#' \item{samechr}{A matrix containing 5 columns as above. Differs from commonmrk in that duplicated markers in either map will have been removed, so all markers map to exactly one chromosome}
#' \item{map1}{The first map - either the originally input object1, or object1$map if it inherits class \code{mpcross}}
#' \item{map2}{The first map - either the originally input object2, or object2$map if it inherits class \code{mpcross}}
#' \item{correlations}{The correlation between positions in map1 and map2 for each chromosome}
#' \item{dup1}{The names of markers duplicated in map1}
#' \item{dup2}{The names of markers duplicated in map2}
#' 
#' Plot produces for a comparison for each chromosome of positions of markers 
#' which are mapped to that chromosome in both maps 
#'
#' Summary function returns printed summary including - number of markers in each map; number of markers common to both maps; number of duplicated markers in each map; number of markers mapped to different chromosomes; correlations between positions on each chromosome.
#' @examples 
#' map1 <- sim.map(len=rep(100, 4), n.mar=51, include.x=FALSE)
#' map2 <- sim.map(len=rep(100, 4), n.mar=51, include.x=FALSE)
#' mc <- mapcomp(map1, map2)
#' summary(mc)
#' plot(mc)

mapcomp <- function(object1, object2) {
	if (!(inherits(object1, "map") | inherits(object1, "mpcross"))) 
		stop("Error: must have map object to work with\n")
	if (!(inherits(object2, "map") | inherits(object2, "mpcross"))) 
		stop("Error: must have map object to work with\n")

	if (class(object1)=="map") map1 <- object1
	if (class(object2)=="map") map2 <- object2
	if (inherits(object1, "mpcross")) map1 <- object1$map
	if (inherits(object2, "mpcross")) map2 <- object2$map

	df1 <- data.frame(mname=unlist(lapply(map1, names)), chr1=rep(names(map1), unlist(lapply(map1, length))), pos1=unlist(map1))
	df2 <- data.frame(mname=unlist(lapply(map2, names)), chr2=rep(names(map2), unlist(lapply(map2, length))), pos2=unlist(map2))

	common <- intersect(df1[,1], df2[,1])

	comm1 <- df1[which(df1[,1] %in% common),]
	comm2 <- df2[which(df2[,1] %in% common),]

	dup1 <- comm1[which(duplicated(comm1[,1])|duplicated(comm1[,1], fromLast=T)),]
	dup1 <- dup1[order(dup1[,1]),]
	dup1[,1] <- as.character(dup1[,1])

	dup2 <- comm2[which(duplicated(comm2[,1])|duplicated(comm2[,1], fromLast=T)),]
	dup2 <- dup2[order(dup2[,1]),]
	dup2[,1] <- as.character(dup2[,1])

	# remove duplicated markers from both maps
	comm1 <- comm1[which(comm1[,1] %in% setdiff(common, c(unique(dup1[,1]), unique(dup2[,1])))),]
	comm2 <- comm2[which(comm2[,1] %in% setdiff(common, c(unique(dup1[,1]), unique(dup2[,1])))),]

	comm1 <- comm1[order(comm1[,1]),]
	comm2 <- comm2[order(comm2[,1]),]

	comm <- cbind(comm1[,1:3], comm2[,2:3])
	
	comm3 <- comm[comm[,2]==comm[,4],]
	comm3 <- comm3[order(comm3[,2], comm3[,3]),]
	comm3[,2] <- as.character(comm3[,2])

	cormaps <- vector(length=length(map1))
	for (i in 1:length(map1)) cormaps[i] <- cor(comm3[comm3[,2]==names(map1)[i], 3], comm3[comm3[,2]==names(map1)[i], 5])

	output <- list()
	output$commonmrk <- comm
	output$samechr <- comm3
	output$map1 <- map1
	output$map2 <- map2
	output$correlations <- cormaps
	output$dup1 <- unique(dup1[,1])
	output$dup2 <- unique(dup2[,1])
	class(output) <- "mapcomp" 
	output
}
#' @S3method plot mapcomp
#' @method plot mapcomp

plot.mapcomp <- function(x, ...)
{
	require(lattice)
	xyplot(pos1~pos2|chr1, data=x$samechr, pch=20, as.table=TRUE, ...)
}

#' @S3method summary mapcomp
#' @method summary mapcomp

summary.mapcomp <- function(object, ...)
{
	cat("Number of markers in map1 is ", length(unlist(object$map1)), "\n")
	cat("Number of markers in map2 is ", length(unlist(object$map2)), "\n")
	cat("Number of common markers is ", nrow(object$commonmrk), "\n")

	cat("Number of duplicated markers in map1 is ", length(object$dup1), "\n")
	cat("Number of duplicated markers in map2 is ", length(object$dup2), "\n")
	cat("Number of markers with differing chromosomes between maps is ", sum(object$commonmrk[,2]!=object$commonmrk[,4]), "\n")
	cat("Correlations between chromosomes are: \n")
	cat(object$correlations, "\n")
}
