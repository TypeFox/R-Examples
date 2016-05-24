library(RUnit)
library(lisp)

results <- runTestSuite(defineTestSuite('lisp',
                                        'lisp',
                                        testFileRegexp='^.+\\.[rR]'))
printTextProtocol(results)
errors <- getErrors(results)
if (errors$nErr || errors$nFail)
  stop('Unit testing failed (vide supra).')
