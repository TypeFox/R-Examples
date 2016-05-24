###
### $Id: softTh.R 22 2014-06-20 20:59:33Z plroebuck $
###

##
## Public
##

##-----------------------------------------------------------------------------
softTh <- function(y, thld) {
    x <- abs(y)
    x <- sign(y) * (x >= thld) * (x - thld)

    x
}

