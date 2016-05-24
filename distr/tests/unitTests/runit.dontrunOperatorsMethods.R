# test for equality of the saved result and the actual result of a dontrun example
test.operatorsMethods <- function() {
  # we compute the actual result
  N <- Norm(0,3)
  P <- Pois(4)
  runit.dontrunOperatorsMethods.actual <- N ^ P
  
  # we load the stored result
  #   we assume that this test is called from within the script in the upper directory
  load("unitTests/runit.dontrunOperatorsMethods.save")
  
  # we compare the stored result with the calculated one
  #   (a comparison with identical (ignoring the environment) gives FALSE...
  result <- all.equal(runit.dontrunOperatorsMethods.actual,
                      runit.dontrunOperatorsMethods.save)
  
  # we check whether the result is TRUE and if not, we write the message
  #   coming from the result
  checkEquals(is.logical(result) && result, TRUE, msg=paste(result, sep="", collapse="\n"))
}

