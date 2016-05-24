"calcrelimp.forboot" <-
function(data,indices,...){
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

## change UG 1.3: account for first column of data including the weights vector
##                and weight the covariance estimation

    dat <- data[indices,-1]
    wt <- data[indices, 1]
    cova <- cov.wt(dat,wt=data[indices,1])$cov
    ausgabe <- list2vec(as(calc.relimp.default.intern(cova,...),"list"))
    return(ausgabe)
}

