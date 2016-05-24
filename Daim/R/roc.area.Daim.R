
roc.area <- function(x, ...) UseMethod("roc.area")



roc.area.Daim <- function(x, method=NULL, col="red", area.color=rgb(1,0,0,alpha=0.5),
		xlab="False positive rate", ylab="True positive rate",
		density=NULL, angle=4, border=NULL, add=FALSE, ...)
{
	if(class(x)[2] != "cv"){
		if(is.null(method))
			method <- "0.632+"
		meth <- charmatch(method,c("0.632+","0.632","loob","sample"))
		all.roc <- FALSE
		if(meth == 4)
			all.roc <- TRUE
		if(!add){
			plot(-1,-1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab,
				main=paste("Method -",method))
		}
		if(all.roc){
			for(i in 1:length(x$all.roc)){
				polygon(c(1,x$all.roc[[i]][,1],0),c(0,x$all.roc[[i]][,2],0),col=area.color, ...)
				lines(x$all.roc[[i]],col=col, ...)
				points(x$all.roc[[i]],col=col, ...)
			}
		}
		if(meth == 1){
			polygon(c(1,1-x$roc$spec632p,0), c(0,x$roc$sens632p,0), col=area.color,
					density=density,angle=angle, border=border, ...)
			lines(1-x$roc$spec632p,x$roc$sens632p,col=col, ...)
			points(1-x$roc$spec632p,x$roc$sens632p,col=col, ...)
		}
		else if(meth == 2){
			polygon(c(1,1-x$roc$spec632,0), c(0,x$roc$sens632,0), col=area.color,
					density=density, angle=angle, border=border, ...)
			lines(1-x$roc$spec632,x$roc$sens632,col=col, ...)
			points(1-x$roc$spec632,x$roc$sens632,col=col, ...)
		}
		else if(meth == 3){
			polygon(c(1,1-x$roc$specloob,0), c(0,x$roc$sensloob,0), col=area.color,
					density=density,angle=angle, border=border, ...)
			lines(1-x$roc$specloob,x$roc$sensloob,col=col, ...)
			points(1-x$roc$specloob,x$roc$sensloob,col=col, ...)
		}
		else if(meth != 4){
			stop("\nThe value of 'method' must be one of '.632+', '.632', 'loob' or 'sample' !\n")
		}
	}
	else{
		if(is.null(method))
			method <- "cv"
		meth <- charmatch(method,c("cv","sample"))
		all.roc <- FALSE
		if(meth == 4)
			all.roc <- TRUE
		if(!add){
			plot(-1,-1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab,
				main=paste("Method -",method))
		}
		if(all.roc){
			for(i in 1:length(x$all.roc)){
				polygon(c(1,x$all.roc[[i]][,1],0),c(0,x$all.roc[[i]][,2],0),col=area.color, ...)
				lines(x$all.roc[[i]],col=col, ...)
				points(x$all.roc[[i]],col=col, ...)
			}
		}
		if(meth == 1){
			polygon(c(1,1-x$roc$specloob,0), c(0,x$roc$sensloob,0), col=area.color,
					density=density,angle=angle, border=border, ...)
			lines(1-x$roc$specloob,x$roc$sensloob,col=col, ...)
			points(1-x$roc$specloob,x$roc$sensloob,col=col, ...)
		}
		else if(meth != 2){
			stop("\nThe value of 'method' must be one of 'cv' or 'sample' !\n")
		}
	}
}



roc.area.Daim.list <- function(x, col="black", area.color=rgb(1,0,0,alpha=0.5),
		xlab="False positive rate", ylab="True positive rate", main="ROC curves",
		density=NULL, angle=4, border=NULL, add=FALSE, ...)
{
	if(!add){
		plot(c(0,1),c(0,1), type="n",xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, main=main)
		grid()
	}
	area.col <- area.color
	line.col <- col
	for(i in 1:length(x)){
		if(length(area.color) > i)
			area.col <- area.color[i]
		if(length(col) > i)
			line.col <- col[i]
		best.cut <- which.max((1-x[[i]]$FPR)+x[[i]]$TPR-1)
		polygon(c(1,x[[i]]$FPR,0),c(0,x[[i]]$TPR,0),col=area.col, 
			density=density,angle=angle, border=border, ...)
		lines(x[[i]]$FPR,x[[i]]$TPR,col=line.col, ...)
		points(x[[i]]$FPR[best.cut], x[[i]]$TPR[best.cut], col=line.col, pch=19, ...)
	}
}


roc.area.Daim.vector <- function(x, col="red", area.color=rgb(1,0,0,alpha=0.5),
		xlab="False positive rate", ylab="True positive rate", main="ROC curve",
		density=NULL, angle=4, border=NULL, add=FALSE, ...)
{
	if(!add){
		plot(c(0,1),c(0,1), type="n",xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, main=main)
		grid()
	}
	best.cut <- which.max((1-x$FPR)+x$TPR-1)
	polygon(c(1,x$FPR,0),c(0,x$TPR,0),col=area.color,
		density=density,angle=angle, border=border, ...)
	lines(x$FPR,x$TPR,col=col, ...)
	points(x$FPR[best.cut], x$TPR[best.cut], col=col, pch=19, ...)
}



roc.area.default <- function(x, ...) {
  stop(paste("Do not know how to handle objects of class", class(x)))
}

