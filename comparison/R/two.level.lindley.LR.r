##########################################################################
# calculates the likelihood ratio for a univariate random effects 
# with between items modelled normal this is taken from Lindley's
# 1977 work and really forms the precursor to the rest of this package
#
#
#
# REQUIRES
#
# control	- a compitem object calculated from the observations from
#		the item considered to be the control item - calculated from
#		two.level.comparison.items() from the file items_two_level.r
# recovered	- a compitem object calculated from the observations from
#		the item considered to be the recovered item - calculated from
#		two.level.comparison.items() from the file items_two_level.r

# background	- a compcovar object calculated from the observavions of the
#		population as a whole - calculated from the two.level.components()
#		function from the file components_two_level.r
#
#
# RETURNS
# 
# LR	- an estimate of the likelihood ratio
##########################################################################
two.level.lindley.LR <- function(control, recovered, background)
{
# first check to make sure that the observations are univariate
if(background@multivariate){stop("Data are multivariate - univariate only allowed for this function")}

# U - is within variance
# C - is between variance
# mu - all means

# redefine some of the items sent
# first the object with the population information
U <- background@v.within
C <- background@v.between
mu <- background@overall.means

# then the objects with the control and recovered
# item information
x <- control@item.means
m <- control@n.replicates
y <- recovered@item.means
n <- recovered@n.replicates


a <- sqrt((1/m) + (1/n))
w <- ((m * x) + (n * y)) / (m + n)

delta.1 <- C + (U / m)
delta.2 <- C + (U / n)
delta.3 <- C + (U / (n + m))
z <- ((delta.2 * x) + (delta.1 * y)) / (delta.1 + delta.2)

bit1 <- (sqrt(delta.1) * sqrt(delta.2)) / (a * sqrt(U) * sqrt(delta.3))
bit2 <- (((x - y)^2) * C) / (a^2 * U * (delta.1 + delta.2))
bit3 <- ((w - mu)^2) / (2 * delta.3)
bit4 <- (((z - mu)^2) * (delta.1 + delta.2)) / (2 * delta.1 * delta.2)
bit5 <- bit4 - bit3

LR <- as.numeric(bit1 * exp(-bit2) * exp(bit5))

return(LR)
}#


