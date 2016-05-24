Mer.Fun <- makeFun(lm(AvgMercury ~ pH, data = FloridaLakes))
Mer.Fun(pH = 7.5) # predicted mercury level at 7.5 pH
Resid <- 1.10 - 0.388 # residual at 7.5 pH
Resid

