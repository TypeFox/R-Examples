## ---- echo=FALSE---------------------------------------------------------
library(memuse)
matmemsize <- function(n) capture.output(memuse::mu(8*2*(n+1)*n))

