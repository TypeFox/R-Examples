#
# Set/Get rounding mode functions for R
# 
# Copyright 2010,2011 Gianluca Amato
#

FE.TONEAREST <- 0
FE.DOWNWARD <- 0x400
FE.UPWARD <- 0x800
FE.TOWARDZERO <- 0xc00

fegetround <- function() {
  .C(r_fegetround, rounding.mode = integer(1))$rounding.mode
}

fesetround <- function(rounding.mode=FE.TONEAREST) {
  .C(r_fesetround, rounding.mode = as.integer(rounding.mode))$rounding.mode
}
