library(bayesPop)
source('test_functions.R')

warn <- options('warn')
options(warn=2)
test.expressions()

# longer tests (comment out for submission)
 # test.prediction()
 # test.prediction.with.prob.migration()
 # test.expressions.with.VE(map=FALSE)
 # test.regional.aggregation()
 # test.life.table()
 
options(warn=warn$warn)