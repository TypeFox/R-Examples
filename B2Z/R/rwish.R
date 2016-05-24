############################################
#This function randomly generates a value  #
#from a Wishart distribution               #
############################################
rwish <-
function (v, S) 
{
#This function is from the package MCMCpack by  
#Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
            (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

