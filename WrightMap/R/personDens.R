personDens <-
function(thetas, yRange = NULL, dim.lab.cex = 0.6, dim.lab.side = 3, dim.lab.adj = 0.5,dim.names = NULL,dim.color = "black",person.points = NULL, person.range = NULL, p.point.col = "black", p.range.col = "gray70",oma = c(0, 5, 0, 5), axis.logits = "Logits",show.axis.logits = TRUE, axis.persons = "Respondents",...) {

	densExt <- function(densElem) {
		xDim <- densElem["y"][[1]]
		yDim <- densElem["x"][[1]]
		xDim <- xDim/max(xDim)

		densInfo <- cbind(xDim, yDim)
		return(densInfo)
	}

	theta.dens <- function(thetas) {
		densList <- apply(thetas, 2, density, na.rm = TRUE)
		distInfo <- lapply(densList, densExt)

		return(distInfo)
	}

	person.plot <- function(distInfo, yRange, xRange, p.points, p.range, p.col, r.col, dim.lab.side, dim.lab.cex, dim.lab.adj, p.cex.lab, p.font.lab, p.lwd) {
		par(mar = c(op$mar[1], 0.2, op$mar[3], 0.1))
		
		plot(distInfo, ylim = yRange, xlim = xRange, type = "l", axes = FALSE, ylab = "", xlab = "", cex.lab = p.cex.lab, font.lab = p.font.lab, lwd = p.lwd, col = attr(distInfo, "dim.color"))
		if(screen() == first)
			mtext(axis.persons, side = 2, line = 1, cex = 0.8, font = 3)
		mtext(attr(distInfo, "dim.name"), side = dim.lab.side, line = -1, cex = dim.lab.cex, font = 1, adj = dim.lab.adj)
		box(bty = "c")
		
		draw.range <- function(upper, lower, col) {
				
				points( c( 0,0), c(lower, upper), type = "l", lwd = 5, lend=2, col = col)
			}
			
		draw.point <- function(pt, col) {
				
				 points(0,pt, pch = 15, cex = .6, col = col)
			}

		 if(!is.null(p.range)) {
				p.range <- matrix(p.range,nrow = 2)
				lower <- p.range[1,]
				upper <- p.range[2,]
				mapply(draw.range,upper,lower,r.col)
			}
			if(!is.null(p.points)) {
				mapply(draw.point,p.points,p.col)
			}
		 
		if (screen() != max(split.screen())) 
			screen(screen() + 1)
	}
	
	thetas <- as.matrix(thetas)
	
	nD <- ncol(thetas)
	if(is.null(yRange))
		yRange <- c(min(thetas),max(thetas))

	xRange <- c(1, 0)
	distInfo <- theta.dens(thetas)
	
	if (is.null(dim.names)) {
		if (!is.null(names(thetas))) {
			dim.names <- names(thetas)
		} else dim.names <- c(paste("Dim", seq(1:nD), sep = ""))
	}
	
	if (ncol(thetas) > 1 & length(dim.color) == 1) {
		dim.color <- rep(dim.color, ncol(thetas))
	}
	
	for (i in 1:nD) {
		attr(distInfo[[i]], "dim.name") <- dim.names[i]
		attr(distInfo[[i]], "dim.color") <- dim.color[i]
	}
	old.screens <- split.screen()
	split.screen(c(1, nD))
	first <- screen()

	op <- par("mar","oma")

	par(oma = oma)

		lapply(distInfo, FUN = person.plot, yRange = yRange, xRange = xRange, p.points = person.points, p.range = person.range, p.col = p.point.col, r.col = p.range.col, dim.lab.cex = dim.lab.cex, dim.lab.side = dim.lab.side, dim.lab.adj = dim.lab.adj, p.cex.lab = 1.3, p.font.lab = 3, p.lwd = 2)
		
		if (show.axis.logits) {
		axis(4, las = 1, cex.axis = 0.7, font.axis = 2)
	}
	mtext(axis.logits, side = 4, line = 1.5, cex = 0.8, font = 3)
	
	curr.screens <- split.screen()
	new.screens <- curr.screens[!(curr.screens %in% old.screens) ]
	close.screen(new.screens)
}
