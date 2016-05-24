###
### $Id: std.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Standard deviation.
###


##-----------------------------------------------------------------------------
std <- function(x, flag = 0) {
    if (flag != 0) {
        stop("biased standard deviation not implemented")
    }

    return(sd(x))
}

