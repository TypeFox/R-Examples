#' Plotting of site-by-species interaction matrices
#'
#' 'Imagine' produces an image plot, grid of small rectangles representing
#' species occurrences in sites, of a given interaction matrix.
#'
#' @param comm community data in the form of a presence absence matrix
#' @param col colors used to plot interactions. First value is the background
#' color (no interaction) and the second color indicates an interaction.
#' @param scores axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores
#' @param order logical. Should the interaction matrix be ordered based on
#' reciprocal averaging scores before plotting?
#' @param fill logical. Should species ranges be made coherent before plotting?
#' @param xlab name of the x axis
#' @param ylab name of the y axis
#' @param yline line that the y-axis label is plotted on.
#' @param xline line that the x-axis label is plotted on.
#' @param sitenames names for each row in the interaction matrix. Default is to
#' not plot names.
#' @param speciesnames names for each site in the interaction matrix. Default
#' is to not plot names.
#' @param binary logical argument indicating whether to ordinate the community
#' matrix based on abundance or binary (default) data.
#' @return Produces an image plot of the interaction matrix. The code is very
#' simple, and may need to be modified if you have long site or species names,
#' or wish to make it prettier than I have the ability to.
#' @author Tad Dallas
#' @export
#' @keywords ordination plotting
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' pres3c=TestMatrices[[6]]
#'
#' #plot interaction matrix
#' Imagine(pres3c, col=c('white','blue'), order=TRUE, fill=FALSE)
#'

Imagine <- function(comm, col=c(0,1), order=TRUE, scores=1, fill=TRUE,  xlab='Species', ylab='Site', yline=2, xline=2, sitenames=rownames(comm), speciesnames=colnames(comm), binary=TRUE){

	require(metacom)
	if(order == TRUE){
	comm <- OrderMatrix(comm, binary= binary, scores=scores)}

	if(fill==TRUE){
	for(i in 1:dim(comm)[2]){
		temp=comm[,i]
		if(sum(temp) < 2){comm[,i]=temp
		}else{
		first <- min(which(temp > 0))
		last <- max(which(temp > 0))
		comm[first:last,i] <- max(temp)
		}
		}
	}

# Format for plotting
 reverse <- nrow(comm) : 1
 comm <- comm[reverse,]

# Image plot
par(mar=c(2,6,6,1))

image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col=col, xlab="", ylab="", axes=FALSE) ; box()


if(length(sitenames)>1){
axis(2, at=1:dim(comm)[1], labels=sitenames, las= 1, cex.axis=1,lwd.ticks=0)
}

if(length(speciesnames)>1){
axis(3, at=1:dim(comm)[2], labels=speciesnames, cex.axis=1, lwd.ticks=0)
}

mtext(xlab, 3, cex=1.5, line=xline)
mtext(ylab, 2, cex=1.5, line=yline)
}
