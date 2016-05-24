#############################################################
#
#	transparentColor <- function(...)
#
#	Internal function allows for defining named colors with opacity
#	alpha = opacity
#
transparentColor <- function(namedColor, alpha = 0.8) {
	res <- col2rgb(namedColor) / 255;
	return(rgb(red = res[1,], green = res[2,], blue = res[3,], alpha = alpha));
}