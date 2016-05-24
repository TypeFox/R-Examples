test.metabolites<-function(){
  # Check output class
  checkTrue(is.vector(metabolites("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(metabolites("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ATP[m_a]","beta-D-Glucose[m_a]","ADP[m_a]","beta-D-Glucose 6-phosphate[m_a]"))
}
