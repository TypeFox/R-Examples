### <======================================================================>
#
# This is a helper function for ecd constructor
# It generates character based classification for an ecd object
#
".ecd.model" <- function(alpha, gamma, sigma, beta, mu, cusp, lambda)
{
    # TODO: Add more classification based on different regions

    
    # standard cusp is defined very strictly
    if( lambda != 3){
        ecd.case.long <- "Special Elliptic"
        ecd.case.short <- "special ecd"
    } else if(alpha == 0 & gamma == 0 & sigma == 1 & beta == 0 & mu == 0){
        ecd.case.long <- "Standard Cusp Elliptic"
        ecd.case.short <- "std csp ecd"
    } else if(cusp > 0){
        ecd.case.long <- "Cusp Elliptic"
        ecd.case.short <- "csp ecd"
    } else if(alpha == 0 & gamma == 0 & sigma != 1 & beta == 0){
        ecd.case.long <- "Cusp Elliptic"
        ecd.case.short <- "csp ecd"
    } else {
        ecd.case.long <- "Elliptic"
        ecd.case.short <- "ecd"
    }
    
    if(beta == 0){
        ecd.case.long.skew <- paste("Symmetric", ecd.case.long, sep = " ")
        ecd.case.short.skew <- paste("symm", ecd.case.short, sep = " ")
    }else{
        ecd.case.long.skew <- paste("Asymmetric", ecd.case.long, sep = " ")
        ecd.case.short.skew <- paste("asymm", ecd.case.short, sep = " ")
    }

    return (unname(c(ecd.case.long.skew, 
                    ecd.case.long,
                    ecd.case.short.skew, 
                    ecd.case.short)))
}
### <---------------------------------------------------------------------->
