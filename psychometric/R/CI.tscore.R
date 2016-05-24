"CI.tscore" <-
function(obs, mx, s, rxx, level=.95)
{
noma <- 1-level
see <- SE.Est(s, rxx)
zs <- - qnorm(noma/2)
mez <- zs*see
that <- Est.true(obs, mx, rxx)
lcl <- that - mez
ucl <- that + mez
mat <- data.frame(SE.Est = see, LCL = lcl, T.Score = that, UCL = ucl)
return(mat)
}

