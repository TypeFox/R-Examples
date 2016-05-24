ncparamF <- function(type1, type2, nu1, nu2)
.C("fpow",  PACKAGE = "fpow", as.double(type1), as.double(type2), as.double(nu1), as.double(nu2), lambda=double(1))$lambda
