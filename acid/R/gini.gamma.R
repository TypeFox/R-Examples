gini.gamma <-
function(p){ # exact from Kleiber p.164
  Gini<-gamma(p+1/2)/(gamma(p+1)*sqrt(pi))
  return(Gini)
}
