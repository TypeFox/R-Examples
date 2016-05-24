## ---- include=FALSE------------------------------------------------------
library(tvm)

## ------------------------------------------------------------------------
rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
rate_curve(fun_r = function(x) rep_len(0.1, length(x)), rate_type = "swap", knots = 1:10)
rate_curve(fun_d = function(x) 1 / (1 + x), knots = 1:10)

## ------------------------------------------------------------------------
r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
r[, c(1, 2)]
r["zero_eff"]
r["swap",c(1.5, 2)]

## ------------------------------------------------------------------------
plot(r)
plot(r, rate_type = "german")
plot(r, rate_type = c("french", "german"))

