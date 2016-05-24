# label an ellipse, allowing a pos argument to specify 
#  center, bottom, left, top, right.  
#  label.pos=NULL uses the correlation to determine top (r>=0) or bottom (r<0) 
#  Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above 
#  and to the right of the max/min coordinates of the ellipse.


label.ellipse <- function(ellipse, label, col="black", 
				label.pos=NULL, xpd=TRUE, 
				tweak=0.5*c(strwidth("M"), strheight("M")), ...){
		
	ellipse <- as.matrix(ellipse)
	if (ncol(ellipse)<2) stop("ellipse must be a 2-column matrix")
#	ranges <- apply(ellipse, 2, range)
	if (is.null(label.pos)) {
		r = cor(ellipse, use="complete.obs")[1,2]
		label.pos <- if (r>0) 3 else 1
	}

	#		index <- if (1:4 %% 2) ... 

	posn <- c("center", "bottom", "left", "top", "right")
	if (is.character(label.pos)) label.pos <- pmatch(label.pos, posn, nomatch=3)-1

	if (label.pos==1) {   # bottom
			index <- which.min(ellipse[,2])
			x <- ellipse[index, 1]
			y <- ellipse[index, 2] + tweak[2]
			}
	else if (label.pos==2) {   # left
			index <- which.min(ellipse[,1])
			x <- ellipse[index, 1] + tweak[1]
			y <- ellipse[index, 2]
			}
	else if (label.pos==3) {   # top
			index <- which.max(ellipse[,2])
			x <- ellipse[index, 1] 
			y <- ellipse[index, 2] - tweak[2]
			}
	else if (label.pos==4) {   # right
			index <- which.max(ellipse[,1])
			x <- ellipse[index, 1] - tweak[1]
			y <- ellipse[index, 2]
			}
	else if (label.pos==0) {   # center
			x <- mean(ellipse[, 1])
			y <- mean(ellipse[, 2]) - tweak[2]
			label.pos <-3
			}
	
	text(x, y, label, pos=label.pos, xpd=xpd, col=col, ...)
}

