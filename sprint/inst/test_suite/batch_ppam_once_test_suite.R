library("RUnit")
library("ff")
library("cluster")
library("golubEsets")


for (nm in list.files("/Users/egrant1/Documents/SPRINT/workspace/sprint/trunk/inst/unitTests/batch/", pattern = "\\.[Rr]$")){
  source(file.path("/Users/egrant1/Documents/SPRINT/workspace/sprint/trunk/inst/unitTests/batch/", nm))
}
library("sprint")
test.suite <- defineTestSuite("batch", dirs = file.path("/Users/egrant1/Documents/SPRINT/workspace/sprint/trunk/inst/unitTests/batch/"),testFileRegexp = 'ppam_once_test.R')

# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

pterminate()
quit()

