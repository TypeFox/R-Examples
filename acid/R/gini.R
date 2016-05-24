gini <-
function(x){
  x.sort<- sort(x)
  n<- length(x)
  ivect<- 1:length(x.sort)
  eins<- rep(1,length(x.sort))
  Gini<- 1/n*(n+1-2*(t(n+1-ivect)%*%x.sort)/(t(eins)%*%x.sort))
  Ginistar <- Gini * n/(n - 1)
  list(Gini=Gini,"bcGini"=Ginistar)
}
