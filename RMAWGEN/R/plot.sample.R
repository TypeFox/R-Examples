NULL 

#' 
#'It makes a plot by sampling (e.g. monthly) the variables \code{x} and \code{y}
#'
#'@param x vector of input data 
#' @param y vector of second input data. Default is \code{normalizeGaussian_severalstations(x=as.data.frame(x),data=as.data.frame(data),origin_x=origin_x,origin_data=origin_data,sample=sample,step=step,prec=prec)[,1]}
#' @param xlim,ylim,xlab,ylab see \code{\link{plot.default}} (Graphic) 
#' @param legend_position  legend position. Default is \code{"topleft"}. See \code{\link{legend}}.
#' @param pch integer single or multi values for \code{pch} (see \code{\link{plot.default}}). Default is 1. 
#' @param col integer single or multi values for \code{col} (see \code{\link{plot.default}}). Default is 1. 
#' @param col_max  maximum value for color scale to apply to \code{\link{rainbow}} or \code{\link{rainbow}}. Utilized if \code{col} is not a vector and both \code{gray} or \code{color} are \code{TRUE}. Default is 0.9 .
#' @param col_min minimum value for color scale to apply to \code{\link{rainbow}} or \code{\link{rainbow}}. Utilized if \code{col} is not a vector and both \code{gray} or \code{color} are \code{TRUE}. Default is 0.1 .
#' @param origin date of the first row of \code{x}. See \code{\link{normalizeGaussian_severalstations}}.
#' @param sample string character containg informatio how to sample \code{x} and \code{y}. Default is NULL. If \code{NULL} no sampling is done.see \code{\link{normalizeGaussian_severalstations}}.  Only \code{NULL} or \code{"monthly"} options are implemented.
#' @param xhist frequency histogram for \code{x}. Default is \code{hist(x,breaks=breaks,plot=FALSE)}. If it is \code{NULL}, no marginal histograms appear.
#' @param yhist frequency histogram for \code{y}. Default is \code{hist(y,breaks=breaks,plot=FALSE)}. If it is \code{NULL}, no marginal histograms appear. =hist(y,breaks=breaks,plot=FALSE),
#' @param axes  see \code{\link{barplot}}
#' @param step,prec see \code{\link{normalizeGaussian_severalstations}}
#' @param breaks see \code{\link{hist}}
#' @param origin_x see \code{\link{normalizeGaussian_severalstations}}. Default value is set equal to \code{origin}.
#' @param origin_data \code{\link{normalizeGaussian_severalstations}}. Default value is set equal to \code{origin}.
#' @param data \code{\link{normalizeGaussian_severalstations}}. Default value is set equal to \code{x}.
#' @param color logical value. If \code{TRUE} and if \code{col} is unspecified, a color scale is applied according to \code{col_min} and \code{col_max} (see \code{\link{rainbow}}). Default is \code{FALSE}.
#' @param gray logical value. If \code{TRUE} and if \code{col} is unspecified, a color scale is applied according to \code{col_min} and \code{col_max} (see \code{\link{gray}}). Default is \code{TRUE}.
#' @param sort logical value. If \code{TRUE}, \code{x} and \code{y} are sorted and a Q-Q plot is presented.  Deafault is \code{FALSE}. 
#' @param valmin_x numerical threshold value over which the variable \code{x} is plotted. It is enabled only if \code{sort} is set \code{TRUE}.
#' @param valmin_y numerical threshold value over which the variable \code{y} is plotted. It is enabled only if \code{sort} is set \code{TRUE}.
#' @param valmin numerical threshold value for \code{valmin_y} and \code{valmin_x} if there are not specified.
#' @param abline arguments for  \code{\link{abline}} function. Default is \code{c(0,1)}. If it is \code{NULL}, \code{\link{abline}} is disabled and not called.
#' 
#' 
#'  @usage
#'  plot_sample(x,
#'  y = normalizeGaussian_severalstations(x = as.data.frame(x), 
#'   data = as.data.frame(data), origin_x = origin_x, origin_data = origin_data, 
#'   sample = sample, step = step, prec = prec)[, 1],
#'  xlim = range(x, na.rm = TRUE),
#'  legend_position = "topleft",
#'  ylim = range(y, na.rm = TRUE), pch = 1, col = 1,
#'  col_max = 0.9, col_min = 0.1, origin, sample = NULL,
#'  xhist = hist(x, breaks = breaks, plot = FALSE),
#'  yhist = hist(y, breaks = breaks, plot = FALSE),
#'  axes = FALSE, step = NULL, prec = 1e-04, breaks = 50,
#'  origin_x = origin, origin_data = origin, data = x,
#'  xlab = "", ylab = "", color = FALSE, gray = TRUE,
#'  sort = FALSE, valmin_x = valmin, valmin_y = valmin,
#'  valmin = -9999, abline = c(0, 1), ...)
#' 
#' 
#' 
#' 
#' @param ...  see graphical parametes on \code{\link{plot.default}} 
#' 
#' @note It makes a plot betwee \code{x} and \code{y} and shows thair respective probibilty histograms. 
#' If \code{y} is missing, it is automatically calculated as one-dimensional Gaussianization of \code{x} through the function \code{\link{normalizeGaussian_severalstations}}.
#' 
#' @return 0 in case of success
#' 
#' @seealso \code{\link{plot.default}},\code{\link{extractmonths}}, see \code{\link{normalizeGaussian_severalstations}}
#' @export
#' @examples 
#' library(RMAWGEN)
#' data(trentino)
#' plot_sample(x=TEMPERATURE_MIN$T0090,sample="monthly",
#'  origin="1958-1-1",axes=FALSE,xlab="Tn [ degC]",
#'  ylab="x")
#' 
#' set.seed(123456) 
#' z <- rexp(10000,rate=0.5) 
#' x <- normalizeGaussian(x=z,data=z) 
#' plot_sample(x=z,xlab="z",ylab="x")
#' 
#' 




plot_sample <- function(x,
		y=normalizeGaussian_severalstations(x=as.data.frame(x),data=as.data.frame(data),origin_x=origin_x,origin_data=origin_data,sample=sample,step=step,prec=prec)[,1],
		xlim=range(x,na.rm=TRUE),
		legend_position="topleft",
		ylim=range(y,na.rm=TRUE),
		pch=1,
		col=1,
		col_max=0.9,
		col_min=0.1,
		origin,
		sample=NULL,
		xhist=hist(x,breaks=breaks,plot=FALSE),
		yhist=hist(y,breaks=breaks,plot=FALSE),
		axes=FALSE,
		step=NULL,
		prec=1e-4,
		breaks=50,
		origin_x=origin,
		origin_data=origin,
		data=x,
		xlab="",
		ylab="",
		color=FALSE,
		gray=TRUE,
		sort=FALSE,
		valmin_x=valmin,
		valmin_y=valmin,
		valmin=-9999,
		abline=c(0,1),
		...) {
	
## CONTROLLARE I VALORI CHE SI PASSANO DA PLOT A PLOT!!!	
	
	if (is.null(xhist)) yhist=NULL
	if (is.null(yhist)) xhist=NULL #  yhist and xhist are always both NULL if one of both is NULL
	
#	if (is.null(y)) {
#		y <- normalizeGaussian_severalstations(x=as.data.frame(x),data=as.data.frame(data),origin_x=origin_x,origin_data=origin_data,sample=sample,step=step,prec=prec)[,1] 
#		xhist=hist(x, breaks=breaks, plot=FALSE)
#		yhist=hist(x, breaks=breaks, plot=FALSE)
#		ylim <- range(y,na.rm=TRUE)
#	}
	
	
	# START BARPLOT 

	if (!is.null(xhist)) {
		def.par <- par(no.readonly = TRUE) # save default, for resetting...
	
	
	
		top <- max(c(xhist$counts, yhist$counts))
		nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(4,1), c(1,4), TRUE)
		layout.show(nf)
		par(mar=c(4,4,1,1))
		plot_sample(x, y, 
				xlim=xlim,
				col=col,
				col_max=col_max,
				col_min=col_min,
				pch=pch,
				ylim=ylim,
				xlab=xlab,
				ylab=ylab,
				sample=sample,
				xhist=NULL,yhist=NULL,
				axes=axes,
				step=step,
				prec=prec,
				breaks=breaks,
				origin_x=origin_x,
				origin_data=origin_data,
				data=data,
				legend_position=legend_position,
				color=color,gray=gray,sort=sort,
				valmin=valmin,
				valmin_x=valmin_x,valmin_y=valmin_y,abline=NULL,...)

		if (!is.null(abline)) abline(abline)
		
		par(mar=c(0,3,1,1))
		barplot(xhist$counts, axes=axes, ylim=c(0,top), space=0)
		par(mar=c(3,0,1,1))
		barplot(yhist$counts, axes=axes, xlim=c(0,top), space=0, horiz=TRUE)
	

		par(def.par)#- reset to default
	
	
	} else if (is.null(sample)) {
	
		xplot <-  x
		yplot <-  y
		if (sort) {
			xplot <- sort(xplot[xplot>valmin_x & !is.na(xplot)]) #/(length(xplot[xplot>valmin_x & !is.na(xplot)])+1)
			yplot <- sort(yplot[yplot>valmin_y & !is.na(yplot)]) #/(length(yplot[yplot>valmin_y & !is.na(yplot)])+1)
		} 
		plot(xplot,yplot,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col,pch=pch,...)
	
	} else if (sample=="monthly") {
		
		months <- months((0.5:11.5)*365/12,abbreviate=TRUE)
		i_months_x <- extractmonths(data=1:length(x),when=months[1],origin=origin_x)
		
		if (length(col)==1) {
			col_min=min(1,max(0,col_min))
			col_max=min(1,max(col_min+0.001,col_max))
			
			if (length(pch)==1) {
				
				
				
				if (color | gray) {
					if (color) col=rainbow(n=length(months),start=col_min,end=col_max)
					if (gray) col=gray((col_max-col_min)*(1:length(months))/(length(months))+col_min)
					pch=array(pch[1],length(months))
				} else {
					
					pch=pch[1]:(pch[1]+length(months)-1)
					col=array(col[1],length(months))

				} 
			} else {
				
				if (color) {
					col=rainbow(n=length(months),start=col_min,end=col_max)
				}
				else if (gray) {
					col=gray((col_max-col_min)*(1:length(months))/(length(months))+col_min)
				} else {
					
					col=array(col[1],length(months))
				
				}
				
			}
		} else if(length(pch)==1) {
		
				pch=array(pch[1],length(months))	
		
		}
			
		
		
		xplot <-  x[i_months_x]
		yplot <-  y[i_months_x]
		if (sort) {
			xplot <- sort(xplot[xplot>valmin_x & !is.na(xplot)]) #/(length(xplot[xplot>valmin_x & !is.na(xplot)])+1)
			yplot <- sort(yplot[yplot>valmin_y & !is.na(yplot)]) #/(length(yplot[yplot>valmin_y & !is.na(yplot)])+1)
		} 			


		plot(xplot,yplot,xlim=xlim,ylim=ylim,col=col[1],pch=pch[1],xlab=xlab,ylab=ylab,...)
		
		for (m in 2:length(months)) {
			
			i_months_x <- extractmonths(data=1:length(x),when=months[m],origin=origin_x)
			xplot <-  x[i_months_x]
			yplot <-  y[i_months_x]
			if (sort) {
				xplot <- sort(xplot[xplot>valmin_x & !is.na(xplot)]) #/(length(xplot[xplot>valmin_x & !is.na(xplot)])+1)
				yplot <- sort(yplot[yplot>valmin_y & !is.na(yplot)]) #/(length(yplot[yplot>valmin_y & !is.na(yplot)])+1)
			} 			
			
			
			
			points(xplot,yplot,pch=pch[m],col=col[m])
			
			
			
		}
		legend(legend_position,pch=pch,col=col,legend=months)
	} else {
		
		xplot <-  x
		yplot <-  y
		if (sort) {
			xplot <- sort(xplot[xplot>valmin_x & !is.na(xplot)]) #/(length(xplot[xplot>valmin_x & !is.na(xplot)])+1)
			yplot <- sort(yplot[yplot>valmin_y & !is.na(yplot)]) #/(length(yplot[yplot>valmin_y & !is.na(yplot)])+1)
		} 
	
	
	
		plot(xplot,yplot,xlim=xlim,ylim=ylim,pch=pch[1],col=col[1],xlab=xlab,ylab=ylab,...)
		
	}
	
	
#	return(0)
	
}


