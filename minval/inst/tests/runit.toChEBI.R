test.to.chEBI<-function(){
  # Check output class
  checkTrue(is.vector(toChEBI("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]")))
  
  # Check output values
  checkEquals(toChEBI("ATP[m_a] + beta-D-Glucose[m_a] => ADP[m_a] + beta-D-Glucose 6-phosphate[m_a]"),"1 15422 + 1 15903 => 1 16761 + 1 17719")
}
