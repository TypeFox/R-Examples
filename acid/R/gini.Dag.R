gini.Dag <-
function(a,p){ # exact from Kleiber p.217
  Gini<-gamma(p)*gamma(2*p+1/a)/(gamma(2*p)*gamma(p+1/a))-1
  return(Gini)
}
