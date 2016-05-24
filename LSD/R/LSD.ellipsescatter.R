

### ellipsescatter ###


#' @export
#' @name ellipsescatter
#' @aliases LSD.ellipsescatter
#' @title Visualize subgroups of two-dimensional data assuming normal distributions
#' @description A scatterplot with additional colored ellipses based on a gaussianity assumption.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param groups a list of indices or vector names to be plotted as groups (not necessarily all of x and y).
#' @param colors a character vector of R build-in colors corresponding to the chosen groups.
#' @param pch the plotting character (to be passed to plot).
#' @param bgcol a R build-in color for non-grouped points.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param xlab x label, standard graphics parameter.
#' @param ylab y label, standard graphics parameter.
#' @param scalesd a numeric value giving the scaling factor for standard deviations in each dimension (defaults to 1).
#' @param level a numeric value (between 0 and 1) giving the confidence level of a pairwise confidence region.
#' @param legend.cex a numerical value giving the amount by which the added legend should be magnified relative to the default.
#' @param location the x and y co-ordinates to be used to position the legend (see 'xy.coords').
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples x = c(rnorm(50),rnorm(100,2),rnorm(50,4))
#' y = (x + rnorm(200,0,0.8))*rep(c(1,4,1),c(50,100,50))
#' x = sign(x)*abs(x)^1.3
#'	
#' groups = list("Green" = 1:50,"Red" = 51:150,"Blue" = 151:200)
#' colors = c("darkgreen","darkred","darkblue")
#' ellipsescatter(x,y,groups,colors,location = "topleft")
#' @keywords scatterplot


ellipsescatter = function(x,y,groups,colors = NULL,pch = 20,bgcol	= "darkgrey",main = "ellipsescatter",xlab = NULL,ylab = NULL,scalesd = 1,level = 0.75,legend.cex = 1,location = "topright",...)
{	
	if(is.null(names(groups))){names(groups) = 1:dim(summary(groups))[1]}
	grouplength = length(names(groups))
	if(is.null(colors)){colors = 1:grouplength}
	ellipse.coordinates = function(r,scale = c(1,1),centre = c(0,0),level = 0.95,t = sqrt(qchisq(level,2)))
	{
		r = min(max(r,-1),1)
		d = acos(r)
		a = seq(0,2*pi,len = 100)
		matrix(c(t*scale[1]*cos(a + d/2) + centre[1],t*scale[2]*cos(a - d/2) + centre[2]),100,2,dimnames = list(NULL,c("x","y")))
	}
	add.ellipses = function(x,y,col)
	{
		points(x,y,col=col,pch=pch)
        points(ellipse.coordinates(round(cor(x,y,method="spearman",use="na.or.complete"),digits=2),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=col,lwd=3)
	}
	if (is.null(xlab)){xlab = deparse(substitute(x))}
	if (is.null(ylab)){ylab = deparse(substitute(y))}
	plot(x,y,pch=pch,col=bgcol,ylab = ylab,xlab = xlab,main="",...)
	mtext(paste(main),3,2,cex=1.25)
	legendtext = c()
	for (k in 1:grouplength){
		add.ellipses(x[groups[[k]]],y[groups[[k]]],col=colors[k])
		legendtext = c(legendtext,paste(names(groups)[k]," (",length(x[groups[[k]]][!is.na(x[groups[[k]]])]),")",sep=""))
	}
	legend(location,legendtext,pt.bg=colors,col=rep("black",grouplength),bg="white",pch=21,cex=legend.cex,inset=0.01)
}	


### aliases ###


LSD.ellipsescatter = ellipsescatter



