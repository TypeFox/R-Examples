#
# Test function predict
#

library(usl)

options(digits = 3)

data(specsdm91)

u <- usl(throughput ~ load, specsdm91)

predict(u)

predict(u, interval = "confidence")
