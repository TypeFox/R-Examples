#
# Test function coef
#

library(usl)

dfr <- data.frame(load=c(1, 2,      4,      6,      8,      10), 
                  tput=c(1, 1.8868, 3.0769, 3.5294, 3.5398, 3.3557))

u <- usl(tput ~ load, dfr)

signif(coef(u)[['sigma']], 3)
signif(coef(u)[['kappa']], 3)
