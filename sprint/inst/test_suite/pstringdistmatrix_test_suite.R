library("RUnit")
library("sprint")
library("Biostrings")
library("stringdist")
library("ff")


for (nm in list.files("../unitTests/pstringdistmatrix/", pattern = "\\.[Rr]$")){
  source(file.path("../unitTests/pstringdistmatrix/", nm))
}

test.suite <- defineTestSuite("pstringdistmatrix", dirs = file.path("../unitTests/pstringdistmatrix/"),testFileRegexp = '*.R')

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
pterminate()
quit()
