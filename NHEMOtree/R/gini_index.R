gini_index <-
function(Klassen_vektor){
  if (length(Klassen_vektor)==0) return(0)
  return(1-sum((table(Klassen_vektor)/length(Klassen_vektor))^2))
}
