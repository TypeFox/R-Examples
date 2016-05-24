library(optimx)
library(svUnit)
optimx.tests <- system.file("unitTests", "runit.all.R", package = "optimx")
cat("Running:", optimx.tests, "\n")
source(optimx.tests)
clearLog()
test.all()
Log()

