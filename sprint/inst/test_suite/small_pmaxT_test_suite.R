library("RUnit")
library("multtest")

for (nm in list.files("../unitTests/pmaxT/", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/pmaxT/", nm))
}
library("sprint")
test.suite <- defineTestSuite("pmaxT", dirs = file.path("../unitTests/pmaxT/"),testFileRegexp = 'runit_small_test.R')

# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

pterminate()
quit()

