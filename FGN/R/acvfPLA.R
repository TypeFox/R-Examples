acvfPLA <-
function(alpha, maxlag){
b <- -1/(2*Reimann(alpha))
c(1,b*(1/(1:maxlag))^alpha)
}
