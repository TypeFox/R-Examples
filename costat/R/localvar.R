localvar <-
function (spec) 
{
#
# Compute time-varying local variance from a spectrum
#
J <- spec$nlevels
n <- 2^J
m <- matrix(spec$D, nrow = J, ncol = n, byrow = TRUE)
#
# Sum spectra over scale
#
answer <- apply(m, 2, sum)
return(answer)
}
