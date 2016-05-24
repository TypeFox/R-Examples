#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[fpScale.R] by DSB Fre 04/07/2008 11:27 (CEST) on daniel@puc.home>
##
## Description:
## Scale an FP term variable appropriately before model search using "fpScale".
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

`fpScale` <-
    function (x, scaling = TRUE) # copy of mfp:::fp.scale
{
    scale <- 1
    shift <- 0
    if (scaling) {
        if (min(x) <= 0) {
            z <- diff(sort(x))
            shift <- min(z[z > 0]) - min(x)
            shift <- ceiling(shift * 10)/10
        }
        range <- mean(x + shift)
        scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
    }
    return(list(shift = shift, scale = scale))
}

