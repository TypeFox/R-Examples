# This function checks that the input values are genuinely integers - 
# see ?is.integer for explanation

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
}


