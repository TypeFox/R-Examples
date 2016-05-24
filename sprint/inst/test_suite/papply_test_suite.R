library("RUnit")
library("ff")

# source the scripts so that the functions are available on all nodes
#for (nm in list.files("sprint-0.2/inst/unitTests/pboot/", pattern = "\\.[RrSsQq]$")) {
#       source(nm)
#}

for (nm in list.files("../unitTests/papply/", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/papply/", nm))
}

test.suite <- defineTestSuite("papply", dirs = file.path("../unitTests/papply/"),testFileRegexp = '*.R')
library("sprint")
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
pterminate()
quit()

