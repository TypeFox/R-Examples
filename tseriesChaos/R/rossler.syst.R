#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
rossler.syst <- function(t, input, parms) {
	x <- input[1]; y <- input[2]; z <- input[3]
	a <- parms[1]; b <- parms[2]; c <- parms[3]

	dx <- -(y+z)
	dy <- x+a*y
	dz <- b+z*(x-c)

	list(c(dx, dy, dz))
}
