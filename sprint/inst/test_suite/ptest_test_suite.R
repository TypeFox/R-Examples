library("RUnit")

for (nm in list.files("../unitTests/ptest/", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/ptest/", nm))
}

test.suite <- defineTestSuite("ptest", dirs = file.path("../unitTests/ptest/"),testFileRegexp = '*.R')
library("sprint")
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
pterminate()
quit()

