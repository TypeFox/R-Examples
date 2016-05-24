#
# Test function summary
#

library(usl)

data(specsdm91)

u <- usl(throughput ~ load, specsdm91)

summary(u, digits=3)
