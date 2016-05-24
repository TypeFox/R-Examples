if(!require("GNE"))stop("this test requires package GNE.")
# ?GNE

xmax <- 10

x <- seq(-xmax, xmax, length=31)


a3Dplot <- function(x, phi, main, ...)
{
	if(!is.function(phi))
		stop("wrong argument phi.")
z <- outer(x, x, phi, ...)
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("red", "orange", "yellow", "green","blue", "blueviolet", "magenta", "pink") ) 
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
par(bg = "white", mar=c(2, 2, 2, 2))
persp(x, x, z, col=color[facetcol], phi=30, theta=65, ticktype="detail",
	main=main, xlab="a", ylab="b", zlab="phi(a,b)")

}

par(mfrow=c(1, 3))
a3Dplot(x, phiFB, main="phiFB, lam->2")
a3Dplot(x, phiKK, lambda=.5, main="phiKK, lam=.5")
a3Dplot(x, phiKK, lambda=.1, main="phiKK, lam=.1")

par(mfrow=c(1, 3))
a3Dplot(x, phiFB, main="phiFB(a,b)")
a3Dplot(x, GrAphiFB, main="Grad a phiFB(a,b)")
a3Dplot(x, GrBphiFB, main="Grad b phiFB(a,b)")

sum( outer(x, x, GrAphiFB) <= 0 )/(31^2)

max(abs(outer(x, x, GrAphiFB)))


par(mfrow=c(1, 3))
a3Dplot(x, phiFB, main="phiFB(a,b)")
a3Dplot(x, phi=phipFB, main="phipFB(a,b), p=1/4", p=1/4)
a3Dplot(x, phi=phipFB, main="phipFB(a,b), p=1", p=1)


par(mfrow=c(1, 3))
a3Dplot(x, GrAphiFB, main="phiFB(a,b)")
a3Dplot(x, phi=GrAphipFB, main="phipFB(a,b), p=1/4", p=1/4)
a3Dplot(x, phi=GrAphipFB, main="phipFB(a,b), p=1", p=1)


