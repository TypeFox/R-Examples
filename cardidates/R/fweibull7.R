`fweibull7` <-
function(x, p) {
   ## add zero for p[7] if missing
   if (length(p) == 6) p <- c(p, 0)
  (p[4] + exp(-(x/p[5])^p[6]))*(1-p[1] * exp(-(x/p[2])^p[3])) - p[7]
}

