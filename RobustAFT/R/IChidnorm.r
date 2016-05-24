IChidnorm <-
function(low,k) {
if (low < -k)            Int <- (pnorm(-k)-pnorm(low)+integrate(Chidnorm,-k,k)$val+1-pnorm(k))/(1-pnorm(low))
if (-k <= low & low < k) Int <- (integrate(Chidnorm,low,k)$val +1- pnorm(k))/(1-pnorm(low))
if (k <= low)            Int <- 1
Int}

