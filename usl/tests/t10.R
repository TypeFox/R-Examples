#
# Test function overhead
#

library(usl)

options(digits = 3, scipen = 6)

data(specsdm91)

overhead(usl(throughput ~ load, specsdm91))

overhead(usl(throughput ~ load, specsdm91),
         newdata = data.frame(load = c(1,2,4)))
