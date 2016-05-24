kernels <- function(xi, x = 0, h = 1, weight=2){
  kweights <- 1/h * switch(weight,
                 "1" =   ifelse(abs(xi-x) <= h, 1-abs(x),0), # Triangular
                 "2" =   ifelse(abs(xi-x) <= h, 3/4*(1- ((xi-x)/h)^2),0), # Epanechnikov
                 "3" =   dnorm((x - xi)/h), # Gaussian
                 "4" =   ifelse(abs(xi-x) <= h, 15/16*(1- ((xi-x)/h)^2)^2,0), # Biweight 
                 "5" =   ifelse(abs(xi-x) <= h, 1/2,0), # Uniform
  )
  kweights
}