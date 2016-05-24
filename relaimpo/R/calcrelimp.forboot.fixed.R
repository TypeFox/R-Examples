"calcrelimp.forboot.fixed" <-
function(data,indices, ...){
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

# yhat and x are be available beforehand
# data is the vector of standardized residuals

    e <- data$e[indices]
    cova <- cov(cbind(data$fit+e,data[,1:(ncol(data)-2)]))

    ausgabe <- list2vec(as(calc.relimp.default.intern(cova,...),"list"))
    return(ausgabe)
}

