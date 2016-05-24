"CI.Rsq" <-
function(rsq, n, k, level=.95)
 {
 noma <- 1-level
sersq <- sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/((n^2-1)*(n+3)))
zs <- - qnorm(noma/2)
mez <- zs*sersq
lcl <- rsq - mez
ucl <- rsq + mez
mat <- data.frame(Rsq = rsq, SErsq = sersq, LCL = lcl, UCL = ucl)
return(mat)
}

