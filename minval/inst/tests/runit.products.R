test.products<-function(){
  # Check output class
  checkTrue(is.vector(products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ADP[m_a]","beta-D-Glucose 6-phosphate[m_a]"))
}
