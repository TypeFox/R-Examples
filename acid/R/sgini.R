sgini <-
function(x,nu=2,w=NULL){
  x   <- sort(x)
  Fx  <- frac.ranks(x,w)
  cov.sam <- cov(x/(mean(x)),(1-Fx)^(nu-1))
  cov.pop <- mean(x/(mean(x))*(1-Fx)^(nu-1))-mean(x/(mean(x)))*mean((1-Fx)^(nu-1))
  Gini     <- -nu*cov.pop 
  Ginistar <- -nu*cov.sam # following Van Kerm sgini Stata function, cov is bias corrected -> bias corrected Gini
  list("Gini"=Gini,"bcGini"=Ginistar)
}
