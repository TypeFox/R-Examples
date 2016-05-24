# test for equality of the saved result and the actual result of a dontrun example
test.minimum <- function() {
  # we compute the actual object
  runit.dontrunMinimum.actual <- Minimum(Norm(), Pois())
  
  # we load the saved object for comparison
  #   we assume that this test is called from within the script in the upper directory
  load("unitTests/runit.dontrunMinimum.save")
  
  # we compare the stored result with the calculated one
  #   (a comparison with identical (ignoring the environment) gives FALSE...
  result <- all.equal(runit.dontrunMinimum.actual,
                      runit.dontrunMinimum.save)
  
  # we check whether the result is TRUE and if not, we write the message
  #   coming from the result
  checkEquals(is.logical(result) && result, TRUE, msg=paste(result, sep="", collapse="\n"))
}

