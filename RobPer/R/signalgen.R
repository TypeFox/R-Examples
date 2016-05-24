signalgen <- function(tt, ytype, pf=1) {
	
    if (ytype=="const"){
        sig <- rep(0,length(tt))
    }
    if (ytype=="sine") {
        sig <- sin(2*pi*tt/pf)
    }
    if (ytype=="trian") {
        ph  <- (tt%%pf)/pf
        sig <- 3*ph*(ph<=2/3)+ 6*(1-ph)*(ph>2/3)
    }
    if (ytype=="peak") {
        ph  <- (tt%%pf)/pf
        sig <- numeric(length(tt))
        sig[ph<=2/3] <- 9*exp(-3*pf^2*(ph[ph<=2/3]-2/3)^2)
        sig[ph>2/3] <- 9*exp(-12*pf^2*(ph[ph>2/3]-2/3)^2)
    }
    return(sig)
}
