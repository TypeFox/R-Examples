# tests for summary stats functions
library(tsna)
require(sna)
require(testthat)

# ---- tests for tSnaStats ----

data(nd_test_nets)
data(moodyContactSim)

test_that('tSnaStats output as expected',{

  # vertex level measure
  output<-tSnaStats(moodyContactSim,'degree')
  
  expect_true(is.ts(output))
  expect_equal(ncol(output),network.size(moodyContactSim))
  expect_equal(nrow(output),756)
  
  # graph level measure
  output<-tSnaStats(moodyContactSim,'gden')
  expect_equal(ncol(output),1)
  expect_equal(nrow(output),756)
  
  # change sampling interval
  output<-tSnaStats(moodyContactSim,'gden',time.interval=100)
  expect_equal(ncol(output),1)
  expect_equal(nrow(output),8)
  
  # check start and end values
  expect_equal(attributes(tSnaStats(moodyContactSim,'gden',start=600,time.interval=50))$tsp,c(600.00, 750.00,   0.02))
  expect_equal(attributes(tSnaStats(moodyContactSim,'gden',end=600,time.interval=200))$tsp,c(40.000, 440.000,   0.005))
  
  expect_error(tSnaStats(nd_test_nets[[1]],'degree'),regexp = "must be a object of class 'networkDynamic'")
  
  # error for non supported function
  expect_error(tSnaStats(moodyContactSim,'foo'),regexp = 'not one of the sna package descriptive statistics currently supported')
  
  # test passing in function args (this would give error if FUN not passed in)
  tSnaStats(moodyContactSim,'centralization',FUN='degree',time.interval = 100)
  
  # test aggregate dur
  dyads<-network.dyadcount(moodyContactSim)
  expect_equal(as.numeric(tSnaStats(moodyContactSim,'gden',start=0,end=400,time.interval = 100,aggregate.dur=100)),c(1/dyads,0,2/dyads,0,2/dyads))
  # this one should miss all the edges
  expect_equal(as.numeric(tSnaStats(moodyContactSim,'gden',start=0,end=400,time.interval = 100,aggregate.dur=0)), c(0,0,0,0,0))
  
})



supported_sna_funs <-c('components',
                       'triad.census',
                       'connectedness',
                       'dyad.census',
                       'efficiency',
                       'gden',
                       'hierarchy',
                       'lubness',
                       'mutuality',
                       # 'centralization', leaving off centralization because it requirs additional args
                       'closeness',
                       'betweenness',
                       #'bonpow',   hidden because has singular value problems
                       'degree',
                       'evcent',
                       'flowbet',
                       'graphcent',
                       'infocent',
                       'loadcent',
                       'prestige')


# chec

test_that('all sna stats used by tSnaStats can be called',{
  
  for (funName in supported_sna_funs){
    message('testing ',funName)
    tSnaStats(moodyContactSim,funName,start=600,time.interval=50)
  }
  
})


# code below used to find which edge cases cause errors
# run tests against all the networkDynamic edge case networks
# reports, but does not trigger errors and warnings (because some of them will give errors)
# for(n in 1:length(nd_test_nets)){
#   for (funName in supported_sna_funs){
#       cat('testing net',n,' (',names(nd_test_nets)[n],') with term ',funName, '\n')
#       # run in tryCacth block so can report all the errors at once
#       tryCatch(
#     {
#       tSnaStats(nd_test_nets[[n]],funName)
#     },
#     warning = function(w){
#       message("\t",w)
#     },
#     error = function(e){
#       message("\t",e)
#     }
#       ) # end try catch
#      cat("\n")
#     }
# }


# ---- tests for tErgmStats ----
require(ergm)

test_that('tSnaStats output as expected',{
  
  output<-tErgmStats(moodyContactSim,'~edges')
  
  expect_true(is.ts(output))
  expect_equal(ncol(output),1)
  expect_equal(nrow(output),756)
  
 
  output<-tErgmStats(moodyContactSim,'~degree(0:3)')
  expect_equal(ncol(output),4)
  expect_equal(nrow(output),756)
  
  # change sampling interval
  output<-tErgmStats(moodyContactSim,'~edges',time.interval=100)
  expect_equal(ncol(output),1)
  expect_equal(nrow(output),8)
  
  # check start and end values
  expect_equal(attributes(tErgmStats(moodyContactSim,'~edges',start=600,time.interval=50))$tsp,c(600.00, 750.00,   0.02))
  expect_equal(attributes(tErgmStats(moodyContactSim,'~edges',end=600,time.interval=200))$tsp,c(40.000, 440.000,   0.005))
  
  expect_error(tErgmStats(nd_test_nets[[1]],'edges'),regexp = "must be a object of class 'networkDynamic'")
  
  expect_error(tErgmStats(as.networkDynamic(network.initialize(3,hyper=TRUE)),'edges'),regexp = "not appropriate for hypergraphic networks")
  
  expect_error(tErgmStats(as.networkDynamic(network.initialize(3,multi=TRUE)),'edges'),regexp = "not appropriate for multiplex networks")
  
  # check that it will add '~' in front if ommited
  tErgmStats(moodyContactSim,'edges',time.interval=200)
  
  # error for non supported function
  expect_error(tErgmStats(moodyContactSim,'foo'),regexp = 'initialization function .+ not found.')
  
  # test aggregate dur
  expect_equal(as.numeric(tErgmStats(moodyContactSim,'edges',time.interval = 100)),c(1,0,2,0,0,0,2,1))
  expect_equal(as.numeric(tErgmStats(moodyContactSim,'edges',time.interval = 100,aggregate.dur=100)),c(1,2,2,0,2,4,8,4))
  
})

# code below used to find which edge cases cause errors
# run tests against all the networkDynamic edge case networks
# reports, but does not trigger errors and warnings (because some of them will give errors)
# does not use a comprehensive set of ergm terms

# testErgmTerms<-c('edges','nsp(1)','balance')
# 
# for(n in 1:length(nd_test_nets)){
#   for (funName in testErgmTerms){
#       cat('testing net',n,' (',names(nd_test_nets)[n],') with term ',funName, '\n')
#       # run in tryCacth block so can report all the errors at once
#       tryCatch(
#     {
#       tErgmStats(nd_test_nets[[n]],funName)
#     },
#     warning = function(w){
#       message("\t",w)
#     },
#     error = function(e){
#       message("\t",e)
#     }
#       ) # end try catch
#      cat("\n")
#     }
# }


