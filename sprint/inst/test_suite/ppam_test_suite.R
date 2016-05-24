library("RUnit")
library("ff")
library("cluster")
library("golubEsets")


for (nm in list.files("../unitTests/ppam/", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/ppam/", nm))
}
library("sprint")
test.suite <- defineTestSuite("ppam", dirs = file.path("../unitTests/ppam/"),testFileRegexp = '*.R')

# === Set up finished ===

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)

pterminate()
quit()

