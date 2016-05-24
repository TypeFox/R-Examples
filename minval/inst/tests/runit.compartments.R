test.compartments<- function(){
  # Check output class
  checkTrue(is.vector(compartments("A[c]")))
  # Check values
  checkEquals(compartments(c("A[c]","B[m]","C[r]")),c("c","m","r"))
}