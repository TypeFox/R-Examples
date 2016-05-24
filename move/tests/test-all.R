library(testthat)
test_check("move")
if(F){
  require(testthat)
  #require(testCoverage)
  tmp<-reportCoverage(sourcefiles=list.files('../R', full.names = T),
                 executionfiles=list.files('testthat', full.names = T))
  tmp<-reportCoverage(sourcefiles=list.files('../R', full.names = T),
                 executionfiles='testthat/test.smallFunctions.R')
  tmp<-reportCoverage(sourcefiles=list.files('../R', full.names = T),
                      executionfiles='test-all.R')
 # commented to prevent test error
  # require(covr);require(shiny)
  cov<-package_coverage(); shine(cov)
}