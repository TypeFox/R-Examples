`fweibull6` <-
function(x, p) {
  (p[4] + exp(-(x/p[5])^p[6]))*(1-p[1] * exp(-(x/p[2])^p[3]))
}

