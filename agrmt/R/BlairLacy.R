BlairLacy <- function(V) {
# measure of concentration, Blair & Lacy 2000, "l"
# argument: V=frequency vector
d <- dsquared(V)^0.5
k <- length(V)
dmax <- ((k-1)/4)^0.5
return(d/dmax)
}
