#'@title Set Colors for a color wheel
#'
#'@description
#'This function set a given number of colors to create a color wheel
#'
#'
#'@param color an R color name or a color in hexadecimal notation
#'@param num integer value indicating how many colors to be added to the wheel
#'@return A character vector with the given color and the set of colors to
#'create a wheel color
#'@author Gaston Sanchez
#'@seealso \code{\link{col2HSV}}
#'@export
#'@examples
#'
#' # create a color wheel based on 'tomato'
#' setColors("tomato", 12)
#' 
#' # set 7 colors for '#3D6DCC'
#' setColors("#3D6DCC", 7)
#'
setColors <-
function(color, num)
{
	# convert to RGB
	rgb_col = col2rgb(color)
	# convert to HSV
	hsv_col = rgb2hsv(rgb_col)[,1]
	# get degree
	hue = hsv_col[1]
	sat = hsv_col[2]
	val = hsv_col[3]
	cols = seq(hue, hue + 1, by=1/num)
	cols = cols[1:num]
	cols[cols > 1] <- cols[cols > 1] - 1
	# get colors with hsv
	colors = hsv(cols, sat, val)
	# transparency
	if (substr(color, 1, 1) == "#" && nchar(color) == 9)
	{
		alpha = substr(color, 8, 9)
		colors = paste(colors, alpha, sep="")
	}
	colors
}
