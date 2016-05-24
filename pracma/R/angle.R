##
##  a n g l e . R
##


Real <- function(z) Re(z)

Imag <- function(z) Im(z)

# Conj <- function(z) Conj(z)

# use abs() for Mod()

angle <- function(z) atan2(Im(z), Re(z))
