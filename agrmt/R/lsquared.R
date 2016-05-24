lsquared <- function(V) {
# measure of concentration, following Blair & Lacy 2000
# argument: V = frequency vector
k <- length(V) # number of categories
lsq <- dsquared(V)/((k-1)/4)
return(lsq)
}
