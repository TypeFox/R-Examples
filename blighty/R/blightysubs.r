# file contains blighty() internal functions
# Copyright - David Lucy January 2006

# function which gets limits which make the plot square
# regardless of the shape of the plot area
# would stand being generalised to producing an arbitary
# aspect ratio for any plot for any plot window
#
# works by calculating new xlimts for any plot window
# to keep even axis within the plot
# it's main use is to get maps and plans with
# the correct aspect ratio - but can be used for
# other reasons
#
# input of two vectors one for the x one for the y
# in the order min x, max x -and- min y, max y
# output is in the same order only the new limits
#
# the output isn't quite accurate I suspect due to margins
# not being accounted for
sqlimits <- function(xlim, ylim)
{

# get the existing x and y limits
x1 <- xlim[1]
x2 <- xlim[2]
y1 <- ylim[1]
y2 <- ylim[2]

frame()
# calculate the existing figure aspect ratio
fig.ratio <- (x2 - x1)/(y2 - y1)
# grab the aspect ratio of the existing plot window
plot.ratio <- par("pin")[1] / par("pin")[2]

# if the x for the plot is larger than that for the
# existing window then we fix the x's and calculate
# new limits for the y
if(fig.ratio >= plot.ratio)
	{
	x1lim <- x1; x2lim <- x2
	ydist <- y2 - y1
	total.ydist <- (ydist * fig.ratio * (1/plot.ratio))
	diff.ydist <- (total.ydist - ydist) / 2
	y1lim <- y1 - diff.ydist; y2lim <- y2 + diff.ydist
	}

# if the x for the plot is smaller than that for the
# existing window then we fix the y's and calculate
# new limits for the x
if(fig.ratio < plot.ratio)
	{
	y1lim <- y1; y2lim <- y2
	xdist <- x2 - x1
	total.xdist <- (xdist * (1/fig.ratio) * plot.ratio)
	diff.xdist <- (total.xdist - xdist) / 2
	x1lim <- x1 - diff.xdist; x2lim <- x2 + diff.xdist
	}

xlims <- c(x1lim, x2lim)
ylims <- c(y1lim, y2lim)

op <- list(xlims, ylims); names(op) <- c("xlims", "ylims")
return(op)
}

