#
# Test when data frame is not ordered
#

library(usl)

options(digits=3)

dfr <- data.frame(load=c(1, 2,      4,      6,      8,      10), 
                  tput=c(1, 1.8868, 3.0769, 3.5294, 3.5398, 3.3557))

dfr <- dfr[order(-dfr[1]), ]

try(u <- usl(tput ~ load, data=dfr))

coef(u)
