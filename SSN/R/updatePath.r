## update path slot in an SpatialStreamNetwork object


updatePath <- function(ssn, filepath){
    ssn@path <- filepath
    print(paste("SSN path updated to", filepath, sep = " "))
    return(ssn)

}
