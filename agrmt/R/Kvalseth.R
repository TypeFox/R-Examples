Kvalseth <- function(V) {
# COV, as cited in Blair & Lacy 2000
# argument: V = frequency vector
COV <- 1-(1-BerryMielke(V))^0.5
return(COV)
}
