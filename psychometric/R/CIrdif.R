"CIrdif" <-
function (r1, r2, n1, n2, level=.95)
 {
rd = r1 - r2
noma <- 1-level
sed <- sqrt((1-r1^2)/n1 + (1-r2^2)/n2)
zs <- - qnorm(noma/2)
mez <- zs*sed
lcl <- rd - mez
ucl <- rd + mez
mat <- data.frame(DifR = rd, SED=sed, LCL = lcl, UCL = ucl)
return(mat)
 }

