#
# Test import of package nlmrt
#

library(usl)

data(specsdm91)

usl(throughput ~ load, specsdm91, method = "nlxb")
