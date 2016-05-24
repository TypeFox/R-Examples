trendHomog <- 
function(x) {

  # functions to extract S and varS
  kens <- function(y) mannKen(y)$S
  kenvars <- function(y) mannKen(y)$varS
  
  # find S and varS for each season
  fr <- frequency(x)
  x1 <-window(x, s = start(x)[1], end = c(end(x)[1],fr), extend = TRUE)
  x1 <- matrix(x1, byrow = TRUE, ncol = fr)
  S <- apply(x1, 2, kens)
  varS <- apply(x1, 2, kenvars)
  
  Z <- S / varS ^ .5
  chi2.tot <- sum(Z ^ 2)
  Zbar <- mean(Z)
  chi2.trend <- fr * Zbar ^ 2
  chi2.homog <- chi2.tot - chi2.trend
  p.value <- pchisq(chi2.homog, fr - 1, 0, FALSE)
  
  list(chi2.trend = chi2.trend, 
       chi2.homog = chi2.homog, 
       p.value = p.value)
}
