# Function to add a simple north-pointer to a blighty() map
# Copyright - David Lucy January 2006

north.pointer <- function(pos="AUTO")
{
# positively associate a value with a missing item
if(missing(pos)){pos <- "AUTO"}

# check for the existence of blighty.mapinfo - stop if not
# the complication comes about because R CMD check generates
# a query if no blighty.mapinfo - so the trick is to create and
# destroy it - then hop out north.pointer() with an appropriate error
if(!exists("blighty.mapinfo"))
{
blighty.mapinfo <- NULL
rm(blighty.mapinfo)
stop("blighty.mapinfo non existent: run blighty() before north.pointer()")
}

# calculate base x and y for automatic selection of position
	if(pos == "AUTO")
		{
		xpos <- ((blighty.mapinfo$xlims[2] - blighty.mapinfo$xlims[1]) * 0.90) + blighty.mapinfo$xlims[1]
		ypos <- ((blighty.mapinfo$ylims[2] - blighty.mapinfo$ylims[1]) * 0.90) + blighty.mapinfo$ylims[1]
		}

# if the user has sent a vector describing where the centre of the "N" should be
	if(is.numeric(pos))
		{
		if(length(pos) != 2){cat("\nWrong type for North pointer coords\n"); stop}
		xpos <- pos[1]
		ypos <- pos[2]
		}

# if the user wishes to visually select the centre of the "N" of the pointer
	if(pos == "select")
		{
		cat("\nUse the pointer to select where the North pointer should be on the map\npress button 1 to select the point, then press button 2 to exit locator\n\n")
		pnts <- locator()
		xpos <- pnts$x
		ypos <- pnts$y
		}

# get a length for the vertical of the north pointer
length <- (blighty.mapinfo$ylims[2] - blighty.mapinfo$ylims[1]) * 0.10
# draw the "N"
text(xpos, ypos, labels="N", adj=0.5)
# Arrow up a north pointer
arrows(xpos, ypos - (1.5 * length), xpos, ypos - (0.2 * length), code=2, angle=20, length=0.1) 
# calculate the length of the cross bar
crossheight <- ypos - (0.85 * length)
crosslength <- length / 15
# draw the cross bar
segments(xpos - crosslength, crossheight, xpos + crosslength, crossheight)
}
