i2qj <- function(i,qmax,jmax){
    if(jmax > qmax)
		stop(paste("Argument \"jmax\" must be less than",
                           "or equal to \"qmax\".\n"))
    M <- matrix(NA,nrow=qmax,ncol=jmax)
    itop <- jmax*(qmax - (jmax-1)/2)
    M[row(M) >= col(M)] <- 1:itop
    q_all <- row(M)[!is.na(M)]
    j_all <- col(M)[!is.na(M)]
    if(missing(i)) {
        i <- 1:itop
    } else if(any(i > itop)) {
            stop("At least one value of \"i\" out of range.\n")
    }
    list(q=q_all[i],j=j_all[i])
}
