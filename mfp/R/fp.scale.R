fp.scale <- function(x, scaling = TRUE)
{
	scale <- 1; shift <- 0
	if(scaling) {
		if(min(x) <= 0) {
			z <- diff(sort(x))
			shift <- min(z[z > 0]) - min(x)
			shift <- ceiling(shift*10)/10
		} 
#
#	range <- median(x+shift)
	range <- mean(x+shift)
#    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
	scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
	}
#
    return(list(shift=shift, scale=scale))
}
