test.orphan.products<-function(){
  # Check output class
  checkTrue(is.vector(orphan.products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(orphan.products("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ADP[m_a]","beta-D-Glucose 6-phosphate[m_a]"))
}
