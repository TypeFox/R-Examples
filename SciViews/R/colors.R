## We often need a red-white-blue color ramp, or a red-yellow-green one
## So, define rwb.colors() and ryg.colors()
rwb.colors <- function (n, alpha = 1, s = 0.9, v = 0.9)
{
	if ((n <- as.integer(n[1L])) <= 0) return(character(0L))
	## Define the initial (red) and final (blue) colors with white in between
	cols <- c(hsv(0, s, v, alpha = alpha),   # Red
			  hsv(0, 0, v, alpha = alpha),   # White
			  hsv(2/3, s, v, alpha = alpha)) # Blue
	## Use a color ramp from red to white to blue
	return(colorRampPalette(cols)(n))
}

## Red-yellow-green palette (take care for color-blind people here)!
ryg.colors <- function (n, alpha = 1, s = 0.9, v = 0.9)
{
	## This is essentially rainbow(), but going from 0 (red) to 2/6 (green)
	return(rainbow(n, s = s, v = v, start = 0, end = 2/6, alpha = alpha))
}

## Slighly different than cm.colors(), allowing for s and v!
## Produce probably better results on a CMYK device (color printer)?
cwm.colors <- function (n, alpha = 1, s = 0.9, v = 0.9)
{
	if ((n <- as.integer(n[1L])) <= 0) return(character(0L))
	## Define the initial (red) and final (blue) colors with white in between
	cols <- c(hsv(1/2, s, v, alpha = alpha), # Cyan
			  hsv(0, 0, v, alpha = alpha),   # White
			  hsv(5/6, s, v, alpha = alpha)) # Magenta
	## Use a color ramp from red to white to blue
	return(colorRampPalette(cols)(n))
}
