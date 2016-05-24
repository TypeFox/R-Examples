IPsidnorm <-
function(low,k) {                                                                  # usata ???
if (low < -k)            Int <- 0
if (-k <= low & low < k) Int <- integrate(Psidnorm,low,k)$val/(1-pnorm(low))
if (k <= low )           Int <- 0
Int}

