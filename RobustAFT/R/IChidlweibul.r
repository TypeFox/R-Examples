IChidlweibul <-
function(r,k) {
c0 <- -0.1351788
if (r < -k+c0)               Int <- (plweibul(-k+c0)-plweibul(r)+integrate(Chidlweibul,-k,k)$val+1-plweibul(k+c0))/(1-plweibul(r))
if (-k+c0 <= r & r < k+c0)   Int <- (integrate(Chidlweibul,r-c0,k)$val +1- plweibul(k+c0))/(1-plweibul(r))
if (k+c0 <= r)               Int <- 1-plweibul(k+c0)
Int}

