### simSummary_unitTests.R
###----------------------------------------------------------------------------

simSummary_unitTests <- structure(
    
  ## --- Function ---
    
function # Run unit tests for package simSummary
  ### Package \pkg{simSummary} has unit tests to test the validity of code behaviour.
  ### Any developeR or useR can easily run these tests as shown in examples bellow. 
() {

  cat("\n\n\n --- Clearing unit test logger --- \n\n\n")
  clearLog()

  cat("\n\n\n --- Running unit tests for package simSummary --- \n\n\n")
  runTest(x=svSuite(tests="package:simSummary"), name="simSummary unit tests")

  ## Possibly run other tests here...

  ## Save log for further processing
  tmp <- Log()

  cat("\n\n\n --- Stats of unit testing (see metadata and summary bellow) --- \n\n\n")
  print(stats(tmp)[, 1:3])

  cat("\n\n\n --- Meta data of unit testing --- \n\n\n")
  print(metadata(tmp))

  cat("\n\n\n --- Summary of unit testing --- \n\n\n")
  print(summary(tmp))

  cat("\n\n\n --- Stats of unit testing (see summary and metadata above) --- \n\n\n")
  print(stats(tmp)[, 1:3])

  ## Make sure we get stop on failure
  cat("\n\n\n")
  errorLog(summarize=FALSE)
      
},## --- Function end ---
    
  ## --- Examples ---

ex=function() {

  ## Run unit tests for package simSummary
  simSummary_unitTests()

})## --- Examples end ---

###----------------------------------------------------------------------------
### simSummary_unitTests.R ends here
