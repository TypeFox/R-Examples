test.chebi.formula<-function(){
  # Check output class
  checkTrue(is.vector(chebi.formula("ATP")))
  # Check output values
  checkEquals(chebi.formula("water"),"H2O")
  # Check case sensitive
  checkEquals(chebi.formula("WATER"),"H2O")
}
