test.is.chebi<-function(){
  # Check output class
  checkTrue(is.vector(is.chebi("ATP")))

  # Check output values
  checkEquals(is.chebi("water"),TRUE)
}