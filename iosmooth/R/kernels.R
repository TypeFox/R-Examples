kernelTrap <- function(x) ifelse(x==0, 1.5/pi, (cos(x) - cos(2*x))/(pi*x^2))
kernelRect <- function(x) ifelse(x==0, 1/pi, sin(x)/(pi*x))