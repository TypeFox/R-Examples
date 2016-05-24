## generating function
NonSymmetric <- function(){ new("NonSymmetric") }

## generating function
EvenSymmetric <- function(SymmCenter = 0){ 
    new("EvenSymmetric", SymmCenter = SymmCenter) 
}

## generating function
OddSymmetric <- function(SymmCenter = 0){ 
    new("OddSymmetric", SymmCenter = SymmCenter) 
}
