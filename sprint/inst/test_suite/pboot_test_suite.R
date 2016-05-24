library("boot")
library("sprint")
library("RUnit")

data(city, package='boot')
data(gravity, package='boot')
data(nuclear, package='boot')
data(aircondit, package='boot')
data(discoveries)
data(trees)
for (nm in list.files("../unitTests/pboot", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/pboot", nm))
}

test.suite <- defineTestSuite("pboot", dirs = file.path("../unitTests/pboot"),testFileRegexp = '*.R')
library("sprint")
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
pterminate()
quit()
