# Function to add a simple bar scale to a blighty() map
# Copyright - David Lucy January 2006

map.scale <- function(pos="AUTO", width=1)
{
# positively associate a value with a missing item
if(missing(pos)){pos <- "AUTO"}

# check for the existence of blighty.mapinfo - stop if not
# the complication comes about because R CMD check generates
# a query if no blighty.mapinfo - so the trick is to create and
# destroy it - then hop out map.scale() with an appropriate error
if(!exists("blighty.mapinfo"))
{
blighty.mapinfo <- NULL
rm(blighty.mapinfo)
stop("blighty.mapinfo non existent: run blighty() before map.scale()")
}



# calculate base x and y for automatic selection of position
	if(pos == "AUTO")
		{
		xpos <- ((blighty.mapinfo$xlims[2] - blighty.mapinfo$xlims[1]) * 0.90) + blighty.mapinfo$xlims[1]
		ypos <- ((blighty.mapinfo$ylims[2] - blighty.mapinfo$ylims[1]) * 0.70) + blighty.mapinfo$ylims[1]
		}

# if the user has sent a vector describing where the centre of the scale should be
	if(is.numeric(pos))
		{
		if(length(pos) != 2){cat("\nWrong type for map scale coords\n"); stop}
		xpos <- pos[1]
		ypos <- pos[2]
		}

# if the user wishes to visually select the centre of the map scale
	if(pos == "select")
		{
		cat("\nUse the pointer to select where the map scale should be on the map\npress button 1 to select the point, then press button 2 to exit locator\n\n")
		pnts <- locator()
		xpos <- pnts$x
		ypos <- pnts$y
		}

# get a length for the scale bar
length <- round((blighty.mapinfo$ylims[2] - blighty.mapinfo$ylims[1]) * 0.10, 0)
vec <- pretty(seq(0, length, length=4))
length <- max(vec)

# draw a simple scale bar with annotation beneath it
segments(xpos - (length / 2), ypos, xpos + (length / 2), ypos, lwd=width)
text(xpos, ypos - (length / 5), labels=paste(length, "km", sep=" "), adj=0.5)
}
