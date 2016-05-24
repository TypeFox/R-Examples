test.unbalanced<-function(){
  # Check output class
  checkTrue(is.vector(unbalanced("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(unbalanced("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),FALSE)
  checkEquals(unbalanced("3 Water[c] => 2 Water[m]"),TRUE)
  checkEquals(unbalanced("Water[r] => beta-D-fructofuranose 1,6-bisphosphate[c_n]"),TRUE)
}
