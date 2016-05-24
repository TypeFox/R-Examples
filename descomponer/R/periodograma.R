periodograma <- function(y) {
# Author: Francisco Parra Rodr?guez
# Some ideas from Gretl 
# http://econometria.wordpress.com/2013/08/21/estimation-of-time-varying-regression-coefficients/ 
cf <- gdf(y)
n <- length(y)
if (n%%2==0) {
m1 <- c(0)
m2 <- c()
for(i in 1:n){
if(i%%2==0) m1 <-c(m1,cf[i]) else m2 <-c(m2,cf[i])}
m2 <-c(m2,0) 
frecuencia <- seq(0:(n/2)) 
frecuencia <- frecuencia-1
omega <- pi*frecuencia/(n/2)
periodos <- n/frecuencia
densidad <- (m1^2+m2^2)/(4*pi)
tabla <- data.frame(omega,frecuencia, periodos,densidad)
tabla$densidad[(n/2+1)] <- 4*tabla$densidad[(n/2+1)]
data.frame(tabla[2:(n/2+1),])}
else {m1 <- c(0)
m2 <- c()
for(i in 1:(n-1)){
if(i%%2==0) m1 <-c(m1,cf[i]) else m2 <-c(m2,cf[i])}
m2 <-c(m2,cf[n]) 
frecuencia <- seq(0:((n-1)/2)) 
frecuencia <- frecuencia-1
omega <- pi*frecuencia/(n/2)
periodos <- n/frecuencia
densidad <- (m1^2+m2^2)/(4*pi)
tabla <- data.frame(omega,frecuencia, periodos,densidad)
data.frame(tabla[2:((n+1)/2),])}
}