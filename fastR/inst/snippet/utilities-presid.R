# partial residual plot
ut.presid <- xyplot(resid(ut.lm1) ~ kwhpday, ut, type=c('p','r'))
