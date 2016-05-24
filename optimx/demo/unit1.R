library(optimx)
library(svUnit)
optimx.test1 <- system.file("unitTests", "runit.1.R", package = "optimx")
cat("Running:", optimx.test1, "\n")
source(optimx.test1)
clearLog()
test.1()
Log()

