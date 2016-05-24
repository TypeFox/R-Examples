gperiodograma <- function(y) {
# Author: Francisco Parra Rodriguez
# Some ideas from Gretl 
# http://econometria.wordpress.com/2013/08/21/estimation-of-time-varying-regression-coefficients/ 
tabla <- periodograma(y)
plot(tabla$frecuencia,tabla$densidad,
main = "Espectro", 
ylab = "densidad",
xlab="frecuencia",type = "l",
col="#ff0000")}