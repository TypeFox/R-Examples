LocalConcavity <- function(psi, dX){

# OUTPUT
# conc : vector with entries conc_k, where conc_k is the difference of slopes
#        to the left and right of psi_k, i.e.
#        conc[k] <= 0 --> psi[k] is a "point of concavity",
#        conc[k] => 0 --> psi[k] is a "point of convexity"
#        where conc_1 = conc_m := 0

m <- length(psi)
conc <- rep(0, m)

# deriv_k := slope between psi_k und psi_(k + 1)
if (m > 2){
    deriv <- (psi[2:m] - psi[1:(m - 1)]) / dX
    conc[2:(m - 1)] <- diff(deriv)
    }
    
return(conc)
}
