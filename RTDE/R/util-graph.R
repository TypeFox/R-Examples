#############################################################################
#   Copyright (c) 2014 Christophe Dutang
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################




simple3Dproj2D <- function(x, y, z, ylim, col="black", posleg="topright", lty=1, ...)
{
    nbz <- length(z)
    
    if(is.matrix(y))
    {
        if(missing(ylim))
            ylim <- range(y[1:nbz,], na.rm=TRUE)
        col <- rep(col, length=nbz)
        lty <- rep(lty, length=nbz)
        
        plot(x, y[1,], type="n", ylim=ylim, lty=lty[1], col=col[1], ...)
    
        for(i in 1:nbz)
            lines(x, y[i,], lty=lty[i], col=col[i])
    }else
    {
        if(missing(ylim))
            ylim <- range(y, na.rm=TRUE)

        plot(x, y, type="l", ylim=ylim, lty=lty[1], col=col[1], ...)
    }
    
	if(is.null(names(z)))
        names(z) <- paste("z=", signif(z,3), sep="")
    
	if(!is.null(posleg))
        legend(posleg, col=col, lty=lty, legend=names(z))
}


qqparetoplot <- function(x, ..., highlight=c("red","cross"))
{
	n <- length(x)
	x <- sort(x)
	highlight <- match.arg(highlight, c("red","cross"))
	
	if(!is.null(names(x)))
        idred <- substr(names(x), 1, 3) == "new"
	else
        idred <- rep(FALSE, n)
    
	p <- (1:n)/(n+1)
	
	plot(-log(1-p[!idred]), log(x[!idred]), ..., xlim=range(-log(1-p)),
        ylim=range(log(x)))
	if(highlight == "red")
        points(-log(1-p[idred]), log(x[idred]), col="red")
	if(highlight == "cross")
        points(-log(1-p[idred]), log(x[idred]), pch="+")
	invisible(list(x=-log(1-p), y=log(x)))
}

