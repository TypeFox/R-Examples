"CI.obs" <-
function (obs, s, rxx, level=.95)
{
noma <- 1-level
sem <- SE.Meas(s, rxx)
zs <- - qnorm(noma/2)
mez <- zs*sem
lcl <- obs - mez
ucl <- obs + mez
mat <- data.frame(SE.Meas = sem, LCL = lcl, OBS = obs, UCL = ucl)
return(mat)
}

