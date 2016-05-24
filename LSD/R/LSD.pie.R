

### LSD.pie ###


#' @export
#' @name LSD.pie
#' @title Custom-built piechart version
#' @description Piecharts at arbitrary position and radii.
#' @param props a numeric vector giving the relations of the pie pieces (need not to be normalized).
#' @param x x-position of the piechart.
#' @param y y-position of the piechart.
#' @param radius a numeric value giving the radius of the piechart.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param edges an integer giving the number of edges the "circle" will have.
#' @param add logical: if \code{TRUE} (\code{FALSE} by default), the pie is added to an existing plot.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param addPercent logical: if \code{TRUE} (\code{FALSE} by default), the percentage of each slice is written inside of the pie.
#' @param textcol a R build-in color for the percentages of \code{addPercent}.
#' @param clockwise if \code{TRUE} (\code{FALSE} by default), slices drawn clockwise (counter clockwise, if \code{FALSE}).
#' @param init.angle a numerical value representing an angle as a starting angle for the drawn slices.
#' @param labels a character vector giving the names for the pie slices.
#' @param cex scaling a numeric value giving the expansion factor for the slice names (if labels are given).
#' @param cex.percentage a numeric value giving the expansion factor for the percentage values (if \code{addPercent = TRUE}).
#' @param border a R build-in color giving the border color (NA by default).
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb, Carina Demel
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples emptyplot(xlim=c(1,9),ylim=c(1,9))
#' mtext(paste("LSD.pie: piecharts"),3,2,cex=1.25)
#' polygon(c(4,2,4,7,8),c(4,8,4,2,8))
#' LSD.pie(sample(1:50,5),4,4,add=TRUE,radius=2,colpal="prgn",alpha=75)
#' LSD.pie(sample(1:50,5),2,8,add=TRUE,colpal="prgn",alpha=75)
#' LSD.pie(sample(1:50,5),8,8,add=TRUE,colpal="prgn",alpha=75)
#' LSD.pie(sample(1:50,5),7,2,add=TRUE,colpal="prgn",alpha=75)
#' @keywords pie


LSD.pie = function(props,x = 0,y = 0,radius = 1,colpal = "prgn",simulate = FALSE,daltonize = FALSE,cvd = "p",edges = 1000,add = FALSE,xlim = c(-1,1),ylim = c(-1,1),main = "LSD.pie: piecharts",alpha = NULL,addPercent = FALSE,textcol = "black",clockwise = FALSE,init.angle=0,labels=c(),cex = 1,cex.percentage = cex,border=NA,...) 				
{
	direction <- if (clockwise) -2 else 2 
	
	# circle coordinates #
	
	xx = sin(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius
	yy = cos(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius
	
	xxtext = sin(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius*0.7
	yytext = cos(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius*0.7
	
	xxlabels = sin(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius*1.25
	yylabels = cos(-(0:edges)/edges*direction*pi+init.angle*pi/180)*radius*1.25
	
	# number of pie pieces #
	
	n = length(props)
	
	# default colors #
	
	cols = colorpalette(colpal,n,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha)
	
	# pie segments #
	
	phi = cumsum(props)/sum(props)*edges
	phi = c(0,round(phi))+1
	
	# plotting region #
	
	if (!add){
		emptyplot(xlim=xlim,ylim=ylim,...)
		mtext(paste(main),3,2,cex=1.25)
	}
	
	# draw pie pieces #
	
	for (j in 1:n){xvec = c(x,x+xx[phi[j]:phi[j+1]])
		yvec = c(y,y+yy[phi[j]:phi[j+1]])
		polygon(xvec,yvec,col="white",border="black")
		polygon(xvec,yvec,col=cols[j],border="black")
	}
	for(j in 1:n){
		if(addPercent){
			xtext = c(x+xxtext[phi[j]:phi[j+1]])
			ytext = c(y+yytext[phi[j]:phi[j+1]])
			text(xtext[(length(xtext)+1)/2], ytext[(length(ytext)+1)/2], paste(round(props[j]/sum(props)*100,1), "%", sep=""), col=textcol, cex=cex.percentage)
		}
		if(length(labels)==n){
			xlabel = c(x+xxlabels[phi[j]:phi[j+1]])
			ylabel = c(y+yylabels[phi[j]:phi[j+1]])
			
			lab = as.character(labels[j])
			if (!is.na(lab) && nzchar(lab)) {
				Px = xlabel[(length(xlabel)+1)/2]
				Py = ylabel[(length(ylabel)+1)/2]
				text(Px,Py,labels[j],col=textcol,cex=cex,adj=ifelse(Px < 0,1,0),xpd=TRUE)
			}
		}
	}
	
	# return invisible #
	
	invisible()
}



