
plot.SCBand <- function(x, y = NULL, xlim = NULL, ylim = NULL, main = NULL, 
	xlab = NULL, ylab = NULL, col = NULL, cex = NULL, pch = NULL, lty = NULL, 
	lwd = NULL, legend = TRUE, where = NULL, text = NULL, legend.cex = NULL, 
	horiz = TRUE, bty = "n", ...) 
{

	object <- x
	xgrid <- seq(min(unlist(object$x)), max(unlist(object$x)), len = object$gridsize)

	if (is.null(xlim)) xlim <- range(unlist(object$x))
	if (is.null(ylim)) {
		ylim <- range(c(unlist(object$y), unlist(y), object$par, object$nonpar,
						object$normscb, object$bootscb), na.rm = TRUE) 
		if (legend) ylim <- extendrange(r = ylim)
	}
	if (is.null(xlab)) xlab <- "Input"
	if (is.null(ylab)) ylab <- "Response"

	ny <- if (is.matrix(y) || is.data.frame(y)) { 1
		} else if (is.list(y) && length(y) == 2) { 2 
		} else if (is.matrix(object$y) || is.data.frame(object$y)) { 1
		} else if (is.list(object$y) && length(object$y) == 2) { 2
		} else 0 
	npar 	<- ifelse(is.null(object$par), 0, 1)
	nnonpar <- if (is.matrix(object$nonpar)) { 2
		} else if (is.numeric(object$nonpar)) { 1
		} else 0 
	nband   <- switch(object$scbtype, no = 0, normal = 1, bootstrap = 1, both = 2)
	nt  	<- ny + npar + nnonpar + nband 

	cex <- if (is.null(cex)) rep.int(.25,ny) else rep_len(cex,ny)
	pch <- if (is.null(pch)) rep.int(1,ny) else rep_len(pch,ny)	
	col <- if (is.null(col)) 1:nt else rep_len(col,nt)
	lty <- if (is.null(lty)) 1:(nt-ny) else rep_len(lty,nt-ny)
	lwd <- if (is.null(lwd)) rep(2,nt-ny) else rep_len(lwd,nt-ny)
	if (is.null(legend.cex)) legend.cex <- 1
	legend.pch <- c(pch,rep(NA,nt-ny)) 
	legend.lwd <- c(rep(1,ny),lwd)
	legend.lty <- c(rep(0,ny),lty)
	col.graph  <- rep(col[(ny+1):nt], rep(c(1,2),c(nt-ny-nband,nband)))
	lty.graph  <- rep(lty, rep(c(1,2),c(nt-ny-nband,nband)))
	lwd.graph  <- rep(lwd, rep(c(1,2),c(nt-ny-nband,nband)))	

	if (is.null(main)) {
		if(npar) {
			main1 <- "Diagnostic for the mean function"
			if(object$model == 0) { 
				main2 <- "Model: no regression effect"
			} else if (object$model == 1) { 
				main2 <- "Model: linear regression" 
			} else if (length(object$model) == 1) { 
				main2 <- paste("Model: polynomial regression of degree",object$model) 
			} else main2 <- paste("Model of dimension", ncol(object$model))
			main <- paste(main1,"\n", main2, sep="")	
		} else if (nnonpar == 2) {
			main = "Comparison of two mean functions"	
		} else {
			main <- "Mean function estimation" }
	}

	if(ny == 0) {
		matplot(xgrid, cbind(object$par, object$nonpar, object$normscb, 
		object$bootscb), type = "l", col = col.graph, lty = lty.graph, lwd = lwd.graph, 
		xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
	} 
	if (ny == 1) {
		matplot(object$x, if(is.null(y)) t(object$y) else t(y), col = col[1],  pch = pch, 
		xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, cex = cex, main = main, ...)			
		matlines(xgrid, cbind(object$par, object$nonpar, object$normscb,
		object$bootscb), col = col.graph, lty = lty.graph, lwd = lwd.graph)			
	}
	if (ny == 2) {
		matplot(object$x[[1]], if(is.null(y)) t(object$y[[1]]) else t(y[[1]]), col = col[1],  
		pch = pch[1], xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, cex = cex[1], 
		main = main, ...)
		matpoints(object$x[[2]], if(is.null(y)) t(object$y[[2]]) else t(y[[2]]), col = col[2],  
		pch = pch[2], cex = cex[2])
		matlines(xgrid, cbind(object$nonpar, object$normscb, object$bootscb), 
		col = col.graph, lty = lty.graph, lwd = lwd.graph)			
	}
	
	if(legend) {
		if (is.null(text))
			text <- c( switch(ny+1, NULL, "Data", c("Data 1","Data 2")), 
			if (ny == 2) { paste("Estimate",1:2) 
			} else if (length(object$model)) { c("Smoothed parametric", "Nonparametric")
			} else "LP estimate", 
			switch(object$scbtype, normal = "Normal SCB", bootstrap = "Bootstrap SCB", 
			both = c("Normal SCB","Bootstrap SCB")))
		if (is.null(where)) where <- "bottomleft"
		legend(where, legend = text, horiz = horiz, pch = legend.pch, cex = legend.cex, 
		col = col, lty = legend.lty, lwd = legend.lwd, bty = bty)
	}		
}

			
	