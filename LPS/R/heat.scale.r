# Draws a heatmap scale for legend
# Author : Sylvain Mareschal <maressyl@gmail.com>
heat.scale <- function(zlim, col.heatmap, at=-10:10, horiz=TRUE, robust=FALSE, customMar=FALSE, title=NA) {
	missing(zlim)
	missing(col.heatmap)
	
	# Background
	if(!isTRUE(customMar)) {
		if(isTRUE(horiz)) { par(mar=c(2.5, 0.5, 2,   0.5))
		} else            { par(mar=c(0.5, 2,   0.5, 4))
		}
	}
	if(isTRUE(horiz)) { plot(x=NA, y=NA, xlim=zlim, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
	} else            { plot(x=NA, y=NA, xlim=0:1, ylim=zlim, xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
	}
	
	# Colors
	breaks <- seq(from=zlim[1], to=zlim[2], along=col.heatmap)
	if(isTRUE(horiz)) { rect(xleft=head(breaks, -1), xright=tail(breaks, -1), ybottom=0, ytop=1, col=col.heatmap, border=NA)
	} else            { rect(xleft=0, xright=1, ybottom=head(breaks, -1), ytop=tail(breaks, -1), col=col.heatmap, border=NA)
	}
	box()
	
	# Unit
	if(isTRUE(robust)) { unit <- "MAD"
	} else             { unit <- "SD"
	}
	
	# Legend
	at <- at[ at != 0 ]
	axis(side=ifelse(isTRUE(horiz), 1, 4), at=at, labels=sprintf("%+i%s", at, ifelse(isTRUE(horiz), "", sprintf(" %s", unit))), las=1)
	axis(side=ifelse(isTRUE(horiz), 1, 4), at=0, labels=ifelse(robust, "median", "mean"), las=1)
	
	# Title
	if(isTRUE(horiz)) {
		if(is.na(title)) title <- sprintf("Gene expression (in %s units)", unit)
		mtext(side=ifelse(isTRUE(horiz), 3, 2), text=title, line=1)
	}
}

