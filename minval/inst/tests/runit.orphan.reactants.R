test.orphan.reactants<-function(){
  # Check output class
  checkTrue(is.vector(orphan.reactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(orphan.reactants("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),c("ATP[m_a]","beta-D-Glucose[m_a]"))
}
