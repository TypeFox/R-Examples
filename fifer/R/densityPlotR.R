##' Given two groups with scores on a quantitative variable, \code{densityPlotR} will draw both distributions, one for each group. 
##'
##' @title Generate a density plot using a formula
##' @param formula a formula object where the grouping variable is on the right side of the ~ and the response variable (quantitative)
##' is on the left side of the ~.
##' @param data a dataset containing the variables listed in the formula 
##' @param colors a vector the same length as the number of levels of the grouping variable indicating the colors to be used
##' for the density lines
##' @export
##' @param ... other arguments passed to the plot function
##' @seealso \code{\link{boxplot}}, \code{\link{prism.plots}}, \code{\link{plotSigBars}}
##' @author Dustin Fife
##' @examples
##' densityPlotR(Petal.Width~Species, data=iris)
densityPlotR = function(formula, data=NULL, colors=NULL, ...){

	#### extract terms
	dv = as.character(formula[[2]])
	iv = as.character(formula[[3]])

	#### get all types
	types = unique(data[,iv])
	
	#### colors
	if (!is.null(colors)){
		if (length(colors)!=length(types)){
			warning("length of colors argument is not the same as the number of groups in the dataset. I will use the default colors.")
			colors = NULL
		}
		colors = string.to.color(as.character(data[,iv]), colors=colors)
	} else {
		colors = string.to.color(as.character(data[,iv]))
	}
	
	#### loop through each type to get range of x and y
	x.em = NA; y.em=NA
	for (i in 1:length(types)){
		subs = which(data[,iv] == types[i])		
		assign(paste("dens", i, sep=""), density(data[subs,dv]))		
		x.em = c(x.em, get(paste("dens", i, sep=""))$x)
		y.em = c(y.em, get(paste("dens", i, sep=""))$y)		
	}
	x.em = na.omit(x.em); y.em=na.omit(y.em)
	
	#### loop through to plot them
	for (i in 1:length(types)){
		subs = which(data[,iv] == types[i])				
		if (i == 1){
			plot(get(paste("dens", i, sep="")), col=colors[subs], xlab=dv, xlim=range(x.em), ylim=range(y.em), ...)
		}
		else {
			lines(get(paste("dens", i, sep="")), col=colors[subs])
		}
	}
	legend("topright", legend=types, text.col=unique(colors), bty="n")
}		
