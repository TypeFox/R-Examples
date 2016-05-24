len.empty.vec <- function(vec, empty) {
# Helper function, computes the number of 'empty' values in a vector vec.
# Intended for remove.genos(). 
# vec: any vector (probably passed by apply() function)
# empty: the symbol appearing in the vector that is considered as an empty value.
#        (could be 0, "0/0", "NA", "1 1", etc)

return(length(which(vec == empty)))

}
