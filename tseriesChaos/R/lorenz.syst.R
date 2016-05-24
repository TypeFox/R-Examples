#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
lorenz.syst <- function(t, input, parms){
	x <- input[1]
	y <- input[2]
	z <- input[3]

	a <- parms[1]
	b <- parms[2]
	c <- parms[3]

	dx <- a * (y - x)
	dy <- x * (b - z) - y
	dz <- x * y + c * z

	list(c(dx, dy, dz))
}
