mmds.2D.plot <- function (x,project = NULL, title = NULL, axis = c(1, 2), xlim = NULL, 
ylim = NULL, outfile.type = NULL,
outfile.name = "mmds",new.plot = TRUE, active.col = x$col[,3], 
active.alpha = 1, sup.col = project$col[,3], active.pch = 20, 
sup.pch = 3, active.lab = FALSE, sup.lab = FALSE, active.cex = 2, 
sup.cex = 2, active.legend.cex = 2, sup.legend.cex = 2, 
active.legend.lwd = 1, sup.legend.lwd = 2, active.lwd = 1, sup.lwd = 4, 
legend = TRUE, active.legend.pos = "bottomleft",
sup.legend.pos = "bottomright", active.legend.name = x$group[,1],
sup.legend.name = project$group[,1], active.legend.col = x$group[,2],
sup.legend.col = project$group[,2], outfile.width = NULL, outfile.height = NULL,
box.lwd = 1, cex.axis = 1, cex.lab = 1, sup.legend.text = 1,
active.legend.text = 1, legend.axis = TRUE, grid = TRUE, axes = TRUE) {

	#check arguments
	if (!inherits(x, "mmds"))
		stop("object of class 'mmds' expected")
	if (any(axis > length(x$eigen.perc)))
		stop("wrong axis")
	if(!is.null(project) && (!inherits(project, "project")))
		stop("object of class 'project' expected")

	#display eigenvalue percentages
	if(legend.axis==TRUE){
	x.lab <- paste("PC", axis[1], " (", round(x$eigen.perc[axis[1]],1), "%)", sep = "")
	y.lab <- paste("PC", axis[2], " (", round(x$eigen.perc[axis[2]],1), "%)", sep = "")
	}
	else{
	x.lab<-""
	y.lab<-""
	}
	if (is.null(xlim)) {
		if (is.null(project)) {
			x.min <- min(x$coord[, axis[1]])
			x.max <- max(x$coord[, axis[1]])
		}
		else {
			x.min <- min(x$coord[, axis[1]], project$coord[, axis[1]])
			x.max <- max(x$coord[, axis[1]], project$coord[, axis[1]])
		}
	xlim <- c(x.min, x.max) * 1.2
	}

	if (is.null(ylim)) {
		if (is.null(project)) {
			y.min <- min(x$coord[, axis[2]])
			y.max <- max(x$coord[, axis[2]])
		}
		else {
			y.min <- min(x$coord[, axis[2]], project$coord[, axis[2]])
			y.max <- max(x$coord[, axis[2]], project$coord[, axis[2]])
		}
	ylim <- c(y.min, y.max) * 1.2
	}
	if(is.null(outfile.type)){
	if (new.plot)
		dev.new()
	}
	else {
	     if(outfile.type=="pdf"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		pdf(paste(outfile.name,"pdf",sep="."),width=outfile.width,height=outfile.height)
		}
	     else if(outfile.type=="png"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		png(paste(outfile.name,"png",sep="."),width=outfile.width,height=outfile.height,units="in",res=150)
		}
	     else if(outfile.type=="tiff"){
		if(is.null(outfile.width)){
			outfile.width<-7
			outfile.height<-7
		}
		tiff(paste(outfile.name,"tiff",sep="."),width=outfile.width,height=outfile.height,units="in",res=150,compression="lzw")
		}
	     else if(outfile.type=="postscript"){
		if(is.null(outfile.width)){
			outfile.width<-0
			outfile.height<-0
		}
		postscript(paste(outfile.name,"ps",sep="."),width=outfile.width,height=outfile.height)
		}
		else {
 		stop("Wrong file format")
}	  
	}
	#print active data
	plot(x$coord[, axis[1]], x$coord[, axis[2]], col = alpha(active.col, active.alpha),
		pch = active.pch, cex = active.cex, lwd = active.lwd, xlab = x.lab, ylab = y.lab, xlim = xlim,
		ylim = ylim,cex.lab=cex.lab,main=title,frame=FALSE,xaxt="n",yaxt="n",mgp=c(1.6,1.6,1.6))
	if (axes == TRUE) {
	  axis(1,lwd=box.lwd,cex.axis=cex.axis)
	  axis(2,lwd=box.lwd,cex.axis=cex.axis)
	}
	if (active.lab)
		text(x$coord[, axis[1]], y = x$coord[, axis[2]],
			labels = rownames(x$coord), pos = 3, col = active.col)

	if (legend)		
		legend (active.legend.pos, bg = "white", active.legend.name, pch = active.pch, pt.cex = active.legend.cex,
		  pt.lwd = active.legend.lwd, col=active.legend.col,title="Active data",cex=active.legend.text,box.lwd=box.lwd)
    if(!is.null(project)) {
    			
		#print supplementary data
		if (legend) 
			legend (sup.legend.pos, bg = "white", sup.legend.name, pch = sup.pch,cex=sup.legend.text,
			  pt.cex = sup.legend.cex, pt.lwd = sup.legend.lwd, col=sup.legend.col,title="Sup data",box.lwd=box.lwd)
		points(project$coord[, axis[1]], project$coord[, axis[2]],
			col = sup.col, pch = sup.pch, cex = sup.cex, lwd = sup.lwd)
		if (sup.lab)
			text(project$coord[, axis[1]], y = project$coord[, axis[2]],
				labels = rownames(project$coord), pos = 3, col = sup.col)
			
	}
	#lastly because in the foreground
	if(grid==TRUE){
	abline (v = 0, lty = 2)
	abline (h = 0, lty = 2)
	}
	box(which="plot",lty="solid",lwd=box.lwd)
	if(!is.null(outfile.type)){
		dev.off()
		print("File created")
	}
}
