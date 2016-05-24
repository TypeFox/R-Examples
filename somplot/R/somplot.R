# somplot: Plotting Kohonen's self-organising maps
# Copyright (C) 2011  Benjamin Schulz and Andreas Dominik
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
# A. Dominik
# TH Mittelhessen
# THM University of Applied Sciences
# Wiesenstrasse 14
# D-35390 Giessen
# andreas.dominik@th-mittelhessen.de
#


som.plot <- function(visfile, datfile, ...) 
{
	coo = read.table(visfile, header=FALSE, fill=TRUE, stringsAsFactors=FALSE, skip = 1, col.names = c("x", "y", "-"))
	dat = read.table(datfile, header=FALSE, fill=TRUE)
	makehexbinplot(data.frame(coo[, c(1,2)], kat = dat[-1, dat[1,1]+1]), ...)
}

hexbinpie <- function(x, y, kat, xbnds=range(x), ybnds=range(y), hbc = NA, pal = NA, hex = "gray", circ = "gray50", cnt = "black", show.counter.border, ...)
{       
	hb  <- hexbin(shape = (diff(ybnds) + 1) / (diff(xbnds) + 1),  x, y, xbnds = xbnds, ybnds = ybnds, IDs = TRUE, xbins = diff(xbnds)*2)
	rx <- 0.5 -> ry
	hexC <- hexcoords(dx = rx, dy = ry / sqrt(3), n = 1)			
	nl <- length(levels(as.factor(kat)))					
	zbnds <- quantile(hb@count, prob = c(0.05, 0.98, show.counter.border), na.rm = TRUE )	# quantile borders for circle diameter and display counter
	zz <- pmax(pmin(sqrt(hb@count / zbnds[2]), 0.85), 0.2)					# circle diameter from 20 to 85% of hexgon diameter
	tt <- unclass(table(kat, hb@cID))							

	for (i in seq(along=zz)) 
	{
		if (!is.na(hex)) 
		{
			polygon(hbc$x[i] + hexC$x, hbc$y[i] + hexC$y, col = NA, border = hex)
		}
    		tp <- pi / 2 - 2 * pi * c(0, cumsum(tt[,i]) / sum(tt[,i]))
		used = FALSE
    		for (j in 1:nl) 
		{
      			if (tp[j+1] == tp[j]) 
			{
				next
			}
			if (j >= 2)
			{
				used = TRUE
	      			pp <- seq(tp[j], tp[j+1], length = floor((tp[j] - tp[j + 1]) * 4) + 2)
	      			xi <- hbc$x[i] + c(0, zz[i] * rx * cos(pp))
	      			yi <- hbc$y[i] + c(0, zz[i] * ry * sin(pp))
	      			polygon(xi, yi, col = pal[j], border = NA, ...)
			}
      		}
		if (!is.na(circ) & used)  
		{
			polygon(hbc$x[i] + rx * zz[i] * cos((1:18) * pi / 9), hbc$y[i] + ry * zz[i] * sin((1:18) * pi / 9), col = NA, border = circ)
		}
    	}
  	for (i in seq(along = zz)) 
	{
    		if ((!is.na(cnt)) & (hb@count[i] > zbnds[3]))
       		{
			text(hbc$x[i], hbc$y[i], hb@count[i], col = cnt, cex = 0.5)
		}
    	}
}

makehexbinplot <-function(data, col = NA, show.legend = TRUE, legend.width = 4, turn = FALSE, window.width = NA, window.height = NA, onlyDefCols = FALSE, scaleX = NA, scaleY = NA, scale = NA, new.xdim = NA, new.ydim = NA, show.box = TRUE, show.axis = FALSE, edit.cols = FALSE, show.counter.border = 0.98, ...)
{
	library(hexbin)

	if (turn)
	{
		data = data.frame(x = data$y, y = data$x, kat= data$kat)
	}

	# scaling of the table
	scaleX = ifelse(!is.na(scale), scale, ifelse(!is.na(scaleX), scaleX, ifelse(!is.na(new.xdim), new.xdim / (diff(range(data$x)) + 1), 1)))
	scaleY = ifelse(!is.na(scale), scale, ifelse(!is.na(scaleY), scaleY, ifelse(!is.na(new.ydim), new.ydim / (diff(range(data$y)) + 1), 1)))

	if (scaleX < 1)
	{
		data$x = floor(data$x * scaleX)
	}
	if (scaleY < 1)
	{
		data$y = floor(data$y * scaleY)		
	}

	if (!show.legend) 
	{
		legend.width = 0
	}

	# calc hbc an fill up empty coordinates with an "empty" element
	pos = 1
	range.x = max(data$x) - min(data$x) + 1
	range.y = max(data$y) - min(data$y) + 1

	hbc = data.frame(x = seq(1,(range.x) * (range.y),1), y = NA)
	for (y in c(min(data$y) : max(data$y)))
	{
		for (x in c(min(data$x):max(data$x)))
		{
			hbc$x[pos] = ifelse(y %% 2 == 1, x + 0.5, x)
			hbc$y[pos] = y * 0.866
			pos = pos + 1
			if (nrow(data[data$x == x & data$y == y,]) == 0)
			{
				data = rbind(data, data.frame(x = x, y = y, kat = ""))
			}
		}
	}

	lvls = levels(as.factor(data$kat))
	lvls = lvls[lvls != ""]

	pal = rainbow(length(lvls))
	if (!is.na(col[1]))
	{
		if (onlyDefCols)
		{
			tmp.pal = rep("white", length(lvls))
		}
		else
		{
			tmp.pal = vector("character", length = length(lvls))
		}
		if (is.data.frame(col))
		{	
			for (i in c(1 : nrow(col)))
			{
				tmp.pal[lvls == col[i,1]] = as.character(col[i,2])
			}
		}
		else
		{
			tmp.pal[c(1:length(col))] = col
		}

		# convert color names into hex values and fill up colors
		if(!onlyDefCols)
		{
			dbl.pal = sprintf("#%02X%02X%02XFF", col2rgb(tmp.pal[tmp.pal != ""])[1,], col2rgb(tmp.pal[tmp.pal != ""])[2,], col2rgb(tmp.pal[tmp.pal != ""])[3,])
			pal = setdiff(pal, dbl.pal)

			for (i in c(1 : length(lvls)))
			{		
				if (is.na(tmp.pal[i]) | tmp.pal[i] == "")
				{
					tmp.pal[i] = pal[1]
					pal = pal[-1]
				}
			}
		}
		pal = tmp.pal
	}
	
	if(edit.cols)
	{
		pal = as.vector(edit(data.frame(kat = lvls, col = pal))[,2])	
	}
	lvls = append("empty", lvls)	
	pal = c("white", pal)	
	
	if(!is.na(window.width))
	{
		window.height = ifelse(is.na(window.height), window.width * (max(hbc$y) - min(hbc$y) - 1 + (range.x / range.y * 2)) / (max(hbc$x) - min(hbc$x) + legend.width), window.height)
		dev.new(width = window.width, height = window.height)
	}
	plot.new()
	plot.window(c(min(hbc$x) - 0.5, max(hbc$x) + 0.5 + legend.width), c(min(hbc$y) - 0.5, max(hbc$y) + 1), asp=0.866)
	if(show.box)
	{
		box()	
	}
	if(show.axis)
	{
		axis(1)
		axis(2)
	}
	if (show.legend)
	{	
		legend("right", lvls[-1], fill=pal[-1], x.intersp=0.2)
	}
	hexbinpie(data$x, data$y, kat=data$kat, hbc = hbc, pal = pal, show.counter.border = show.counter.border, ...) 
}
