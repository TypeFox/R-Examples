#'@title col2HSV: converts a color to HSV in hexadecimal notation
#'
#'@description
#'col2HSV converts an R color (or a set of colors) into an HSV color model, and
#'then returns the color names in hexadeciaml notation
#'
#'@param color an R color name or a color in hexadecimal notation
#'@return A character vector with the color(s) name(s) in hexadecimal notation
#'@author Gaston Sanchez
#'@seealso \code{\link{wheel}}
#'@export
#'@examples
#'
#' # convert 'tomato'
#' col2HSV("tomato")
#'
col2HSV <-
function(color)
{
	# convert to RGB
	rgb_col = col2rgb(color)
	# convert to HSV
	hsv_col = rgb2hsv(rgb_col)
	if (length(color) == 1)
	{
		# get degree
		hue = hsv_col[1]
		sat = hsv_col[2]
		val = hsv_col[3]
		# get colors with hsv
		hex_col = hsv(hue, sat, val)
	}
	if (length(color) > 1)
	{
		hex_col = rep("", length(color))
		for (j in 1:length(color))
		{
			hex_col[j] = hsv(hsv_col[1,j], hsv_col[2,j], hsv_col[3,j])
		}
	}
	hex_col
}
