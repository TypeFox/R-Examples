## Plotting probability distribution

plotNodeProb <- function(x0, y0, x1, y1, prob, seglist, state, cpal, pruned, index, horiz=TRUE, 
	axes=c("no", "no"), bgcol="grey95", pruned.col="red", cex.axes=0.6, by.state=FALSE, type="b", 
	lwd=par("lwd"), x.prop=FALSE, frame=TRUE, ygrid=by.state, cm.pal="blue", xsrt=0, ysrt=0) {

	if (getOption("verbose")) {
		cat("    [-] plotNodeProb: x0=", x0, ", y0=", y0, ", x1=", x1, ", y1=", y1, "\n")
	}

	A <- colnames(prob)

	if (!horiz) {
		x0.tmp <- x0
		x0 <- y0
		y0 <- x0.tmp
		x1.tmp <- x1
		x1 <- y1
		y1 <- x1.tmp
	}	

	xsize <- x1-x0
	ysize <- abs(y1-y0)

	nbseg <- nrow(seglist)
	seg.lab <- rownames(seglist)
	has.segment <- rownames(seglist) %in% rownames(index)	

	## Filling the plotting area with the background color
	rect(x0, y0, x1, y1, col=bgcol)

	pdim <- xsize/nbseg
	xleft <- x0

	if (type=="a") { 
		by.state <- TRUE 
		ygrid <- TRUE
	} else if (type=="cm") {
		by.state <- TRUE 
		ygrid <- FALSE

		if (cm.pal=="blue") {
			## cm.map <- rev(sequential_hcl(10, h = 210, c = c(80, 30), l = c(20, 80), power = 1.5))
			cm.map <- c("#9BCFDA", "#91C9D6", "#7FBFCD", "#65B2C1", "#40A3B4", "#0093A5", "#008195",
 				"#006F85", "#005D75", "#004E69")
		} else if (cm.pal=="heat") {
			cm.map <- c("#E2E6BD", "#E5DB80","#E4C96A","#E0B458","#D89E4B","#CE8642","#C16D3E",
				"#B1523C", "#A0353C", "#8E063B")
		} else if (cm.pal=="heat20") { 
			cm.map <- c("#E2E6BD","#E4E38E","#E5DC81","#E5D476","#E4CB6C","#E3C263","#E1B85B",
				"#DEAD53","#DAA24D","#D69748","#D18C44","#CB8141","#C5753F","#BE693E","#B75C3D",
				"#B04F3C","#A8423C","#9F333C","#97223C","#8E063B")
		}
	}

	if (by.state) {
		## stsep <- (0.10*ysize)/length(A)
		stsep <- 0
		stsize <- (ysize-((length(A)-1)*stsep))/length(A)
		ytmp <- y0
	} else { 
		stsize <- ysize
	}

	if (ygrid) {	
		for (s in 1:(length(A)-1)) {
			ytmp <- ytmp-stsize
			segments(x0, ytmp, x1, ytmp, col="grey50")
		}
	}

	## 
	for (g in 1:nbseg) {
		xright <- xleft+pdim
		ytmp <- y0

		if (has.segment[g]) {
			idseg <- rownames(seglist)[g]

			for (s in 1:length(A)) {
				p <- prob[idseg, s]
				ybot <- ytmp-(p*stsize)

				if (type=="b") {
					rect(xleft, ytmp, xright, ybot, col=cpal[s], border=NA)
					ytmp <- if (by.state) { ytmp-stsize } else { ybot }
				} else if (type=="s") {
					segments(xleft, ybot, xright, ybot, col=cpal[s], lwd=lwd)

					if (g>1) { 
						yprev <- ytmp-(prob[rownames(seglist)[g-1], s]*stsize)
						segments(xleft, yprev, xleft, ybot, col=cpal[s], lwd=lwd) 
					}
					ytmp <- if (by.state) { ytmp-stsize } else { y0 }
				} else if (type=="l") {
					if (g>1) { 
						yprev <- ytmp-(prob[rownames(seglist)[g-1], s]*stsize)
						segments(xleft-(pdim/2), yprev, xleft+(pdim/2), ybot, col=cpal[s], lwd=lwd) 
					}
					ytmp <- if (by.state) { ytmp-stsize } else { y0 }
				} else if (type=="a") {
					xc <- xleft+(pdim/2) 
					xd <- if (x.prop) { (pdim*p)/2 } else { pdim/2 }
					yc <- ytmp-(stsize/2)
					yd <- (stsize*p)/2
					rect(xc-xd, yc+yd, xc+xd, yc-yd, col=cpal[s], border=NA)
					ytmp <- ytmp-stsize
				} else if (type=="cm") {
					rect(xleft, ytmp, xright, ytmp-stsize, col=cm.map[floor(p*10)+1], border=NA)
					ytmp <- ytmp-stsize
				}
			}
		}
		xleft <- xright
	}	

	## Frame surrounding the area
	if (frame) { rect(x0, y0, x1, y1) }

	## Plotting the axes
	## x axis
	if (axes[1]=="bottom") {
		axe.offset <- if (nbseg>1 && any(pruned, na.rm=TRUE)) { ysize*0.3 } else {0}
		segments(x0+(pdim/2), y0, x0+(pdim/2), y0+(0.10*ysize))
		segments(x1-(pdim/2), y0, x1-(pdim/2), y0+(0.10*ysize))
		text(x=c(x0+(pdim/2),x1-(pdim/2)), y=y0+(0.25*ysize), labels=c(1,nbseg), cex=cex.axes, srt=xsrt)
	} else if (axes[1]=="top") {
		segments(x0+(pdim/2), y1, x0+(pdim/2), y1-(0.10*ysize))
		segments(x1-(pdim/2), y1, x1-(pdim/2), y1-(0.10*ysize))
		text(x=c(x0+(pdim/2),x1-(pdim/2)), y=y1-(0.25*ysize), labels=c(1, nbseg), cex=cex.axes, srt=xsrt)
	}

	## y axis
	if (!by.state) {
		if (axes[2]=="left") {
			## segments(x0-(0.1*xsize), y0, x0-(0.1*xsize), y1)
			segments(x0, y0, x0-(0.10*xsize), y0)
			segments(x0, y1, x0-(0.10*xsize), y1)
			segments(x0, (y0+y1)/2, x0-(0.05*xsize), (y0+y1)/2)
			text(x=c(x0-(0.25*xsize),x0-(0.25*xsize)), y=c(y0, y1), labels=c(0,1), cex=cex.axes, srt=ysrt)
		} else if (axes[2]=="right") {
			segments(x1+(0.1*xsize), y0, x1+(0.1*xsize), y1)
			segments(x1+(0.1*xsize), y0, x1+(0.15*xsize), y0)
			segments(x1+(0.1*xsize), y1, x1+(0.15*xsize), y1)
			text(x=c(x1+(0.25*xsize),x1+(0.25*xsize)), y=c(y0, y1), labels=c(0,1), cex=cex.axes, srt=ysrt)
		}
	}

	## A bar showing the pruned and unpruned nodes
	if (nbseg>1 && any(pruned, na.rm=TRUE)) {
		rect(x0, y0+(ysize*0.1), x1, y0+(ysize*0.3), col=bgcol, border=NA)
		xleft <- x0
		for (g in 1:nbseg) {
			xright <- xleft+pdim
			if (has.segment[g] %in% rownames(prob)) {
				if (pruned[has.segment[g]]) {
					rect(xleft, y0+(ysize*0.1), xright, y0+(ysize*0.3),
						col=pruned.col, border=NA)
				} else if (!pruned[has.segment[g]]) {
					rect(xleft, y0+(ysize*0.1), xright, y0+(ysize*0.3), 
						col="green", border=NA)
				}
			}
			xleft <- xright
		}
		rect(x0, y0+(ysize*0.1), x1, y0+(ysize*0.3))
	} else if (nbseg==1 && pruned) {
		segments(x0, y0, x1, y1, col = pruned.col)
	}
}


 
