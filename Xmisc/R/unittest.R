
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Aug 06 14:57:47 EDT 2014 -0400 (Week 31)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


##' Unit testing for developing R packages
##'
##' 
##' @title Unit testing for developing R packages
##' @description
##' Unit testing for developing R packages
##' @author Xiaobei Zhao
##' @field pkg character the name of the package
##' @field testDpath character the absolute directory names where to
##' look for test files. Default: <pkg>/tests
##' @field testFnameRegexp character Regular expression for matching
##' test file names. Default: *.R
##' @field testFuncRegexp character Regular expression for matching
##' test functions. Default: test.*
##' @examples
##' \dontrun{
##' pkg <- 'Xmisc'
##' test.obj <- UnitTest$new(pkg=pkg)
##' test.obj$runme()
##' }
##' @exportClass UnitTest
UnitTest <- 
  setRefClass(
    'UnitTest',
    list(
      pkg='character',
      testDpath='character',
      testFnameRegexp='character',
      testFuncRegexp='character',
      out='list'
      ),
    contains='xRefClass',
    methods=list(
      initialize=function(...){
        .idx <- c(pkg=1)
        callSuper(...,.index=.idx)
        setme()
      },
      setme=function(){
        if (!length(testDpath)){
          testDpath <<- system.file("tests", package=pkg)
        }
        if (!is.dir(testDpath)){
          stop('UnitTest | testDpath is not available. (', testDpath, ')')
        }
        
        if (!length(testFnameRegexp)){
          testFnameRegexp <<- sprintf('^.+\\.R$')
        }
        if (!length(testFuncRegexp)){
          testFuncRegexp <<- sprintf('^test.+')
        }
      },
      defineme=function(){
        ## check availability of RUnit
        if (!require("RUnit", quietly=TRUE)) {
          warning(pkg,' | R package `RUnit` must be available for unit test.')
          return()
        } else {
          library(package=pkg, character.only=TRUE)
          test.suite <- RUnit::defineTestSuite(
            name=pkg,
            dirs=file.path(testDpath),
            testFileRegexp=testFnameRegexp,
            testFuncRegexp=testFuncRegexp
            )
        }
        out$test.suite <<- test.suite
      },
      runme=function(){
        if (!length(out$test.suite)){
          defineme()
        }
        if (!length(out$test.suite)){
          return()
        }
        
        test.result <- RUnit::runTestSuite(out$test.suite)
        out$test.result <<- test.result
      },
      printme=function(){
        if (!length(out$test.result)){
          runme()
        }
        if (!length(out$test.result)){
          return()
        }
                
        RUnit::printTextProtocol(out$test.result)
        tmp <- RUnit::getErrors(out$test.result)
        if(tmp$nFail > 0 | tmp$nErr > 0) {
          stop(
            pkg,' | Unit testing failed (',
            tmp$nFail, ' test failure(s); ', tmp$nErr,' R error(s).'
            )
        }
        invisible()
      }
      )
    )
