`getdata` <-
function(imwd, switch, co.type = "sqr", verbose = FALSE)
{
#
# This routine extracts the detail coefficients from an imwd(tpe="station") object
# and instead of storing them in a wLNx format, stores them in a 3-D array.
# 
# Reason for doing this: I find it much easier to think of accessing an array.
#
#
#
#       Check class of imwd to amke sure that we can compute the DDEWS
#
now <- proc.time()[1.:2.]
ctmp <- class(imwd)
if(is.null(ctmp))
   stop("imwd has no class")
else if(ctmp != "imwd")
   stop("imwd is not of class imwd")
#
# Are we still going ... ok, then let us proceed:
#
if(imwd$type == "wavelet") 
   stop("Decomposition  not performed using the Non-Decimated Wavelet Transfrom")
Csize <- 2.^(imwd$nlevels)
first.last.d <- imwd$fl.dbase$first.last.d
first.last.c <- imwd$fl.dbase$first.last.c
J <- imwd$nlevels
firstD <- first.last.d[J, 1.]
lastD <- first.last.d[J, 2.]
LengthD <- lastD - firstD + 1.
if(switch == "level") {
   unddews <- array(0., dim = c(3. * J, LengthD, LengthD))
   for(level in (J:1.)) {
      ndata <- 2.^(level - 1.)
      firstD <- first.last.d[level, 1.]
      lastD <- first.last.d[level, 2.]
      LengthD <- lastD - firstD + 1.
#####################################################
# Extract CD for this level
#####################################################
      sel <- seq(from = (1. - firstD), length = ndata)
      nm <- lt.to.name(level - 1., "CD")
#########################################################
# Extract DC for this level
######################################################
      msub1 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
      nm <- lt.to.name(level - 1., "DC")
#
############################################################
# Extract DD for this level
#########################################################
      msub2 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
      nm <- lt.to.name(level - 1., "DD")
#####
#
# However, as Herrick has noted (personal correspondence)
# DC is actually CD whereas CD is DC.
#
# Therefore CD is actually the horizontal component
# and similarly DC is the vertical component.
# So msub1 contains the horizontal coefficients, 
# whilst msub2 contains the vertical coefficients.
#####
#############################################################
#       Work out if we want to display the absolute values or the actual
#       values. By default, the coefficients we extract are the square of the coeffs
# returned by the NDWT (these are the ones used to calculate the CDDEWS).
#However, being nice, options are included so that one can see the abs value, - abs value
#and the raw coefficient values.
##########################################################
      msub3 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
# Vertical
      if(co.type == "sqr") {
         msub1 <- msub1^2.
         msub2 <- msub2^2.
         msub3 <- msub3^2.
      }
#
# Horiztonal
      unddews[(3. * (J - level) + 1.),  ,  ] <- msub2
#
# Diagonal
      unddews[(3. * (J - level) + 2.),  ,  ] <- msub1
      unddews[(3. * (J - level) + 3.),  ,  ] <- msub3
   }
}
if(switch == "direction") {
   unddews <- array(0., dim = c(3. * J, LengthD, LengthD))
   for(level in (J:1.)) {
      ndata <- 2.^(level - 1.)
      firstD <- first.last.d[level, 1.]
      lastD <- first.last.d[level, 2.]
      LengthD <- lastD - firstD + 1.
#####################################################
# Extract CD for this level
#####################################################
      sel <-seq(from = (1. - firstD), length = ndata) 
      nm <- lt.to.name(level - 1., "CD")
#
#########################################################
# Extract DC for this level
######################################################
      msub1 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
      nm <- lt.to.name(level - 1., "DC")
#
############################################################
# Extract DD for this level
#########################################################
      msub2 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
      nm <- lt.to.name(level - 1., "DD")
#####
#
# Therefore CD is actually the horizontal component
# and similarly DC is the vertical component.
# So msub1 contains the horizontal coefficients, 
# whilst msub2 contains the vertical coefficients.
#####
      msub3 <- matrix(imwd[[nm]], nrow = LengthD, ncol = LengthD)
# Vertical
      if(co.type == "sqr") {
         msub1 <- msub1^2. 
         msub2 <- msub2^2.
         msub3 <- msub3^2.
      }
# Horizontal
      unddews[(J - level + 1.),  ,  ] <- msub2
# Diagonal
      unddews[(2. * J - level + 1.),  ,  ] <- msub1
      unddews[(3. * J - level + 1.),  ,  ] <- msub3
   }
}
unddews
}

