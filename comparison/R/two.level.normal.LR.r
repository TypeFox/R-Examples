##########################################################################
# calculates the likelihood ratio for a mult-variate random effects 
# with between items modelled normal rather than kernal
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
two.level.normal.LR <- function(control, recovered, background)
{

# redefine some of the items sent
# first the object with the population information
U <- background@v.within
C <- background@v.between
mu <- background@overall.means

# then the objects with the contro and recovered
# item information
cont.means <- control@item.means
n.cont <- control@n.replicates
rec.means <- recovered@item.means
n.rec <- recovered@n.replicates

# weighted mean for the control and recovered item
y.star <- ((n.cont * cont.means) + (n.rec * rec.means)) / (n.cont + n.rec)

# calculate some of the repeated components
diff.cont.rec <- cont.means - rec.means
diff.cont.mu <- cont.means - mu
diff.rec.mu <- rec.means - mu
diff.y.star.mu <- y.star - mu


# This code independently written by Agnieszka Rzepecka of the IFR, and tidied up a bit by David Lucy
# output agrees with previous efforts - thus we develop some confidence it's doing the right thing
# I keep on expecting it to fall over in the univariate case but it seems not to so we'll use it
# until it starts to be a problem
# Numerator calculation
nom1 <- exp(-1/2*t(diff.cont.rec) %*% solve(U/n.cont+U/n.rec) %*% (diff.cont.rec))*(det(U/n.cont+U/n.rec))^(-1/2)
nom2 <- exp(-1/2*t(diff.y.star.mu) %*% solve(U/(n.cont+n.rec)+C) %*% (diff.y.star.mu))*(det(U/(n.cont+n.rec)+C))^(-1/2)
nom <- nom1 * nom2

# Denominator calculation
denom1 <- exp(-1/2*t(diff.cont.mu) %*% solve(U/n.cont+C) %*% (diff.cont.mu))*(det(U/n.cont+C))^(-1/2)
denom2 <- exp(-1/2*t(diff.rec.mu) %*% solve(U/n.rec+C) %*% (diff.rec.mu))*(det(U/n.rec+C))^(-1/2)
denom <- denom1 * denom2

LR <- as.numeric(nom / denom)

#####################################################################################################
# original code from 2004 which follows the Aitken & Lucy paper (p.115) closely #####################
# this is here so people trying to follow the code with the paper can do so #########################
#####################################################################################################
#D1 <- U / n.cont
#D2 <- U / n.rec
#H2 <- t(y.star - mu) %*%  solve(U/(n.rec + n.cont) + C) %*% (y.star - mu)
#H3 <- t(diff.cont.rec) %*% solve(D1 + D2) %*% (diff.cont.rec)
#num1 <- exp(- 0.5 * (H2 + H3))
#num2 <- sqrt(det(2 * pi * solve(((n.cont + n.rec) * solve(U) ) + solve(C))))
#num <- num1 * num2
#mu.star <- solve(solve(D1 + C) + solve(D2 + C)) %*% ((solve(D1 + C) %*% cont.means) + (solve(D2 + C) %*% rec.means))
#H4 <- t(mu - mu.star) %*% (solve(D1 + C) + solve(D2 + C)) %*% (mu - mu.star)
#H5 <- t(diff.cont.rec) %*% solve(D1 + D2 + (2 * C)) %*% diff.cont.rec
#den1 <- 1 / sqrt(det(2 * pi * C))
#den2 <- sqrt(det(2 * pi * solve(n.cont * solve(U) + solve(C))))
#den3 <- sqrt(det(2 * pi * solve(n.rec * solve(U) + solve(C))))
#den4 <- exp(-0.5 * (H4 + H5))
#lr <- num / prod(den1, den2, den3, den4)
######################################################################################################

return(LR)
}#

