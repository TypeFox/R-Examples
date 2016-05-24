NAtoZero <- function(v, value = 0){
    v[is.na(v) == TRUE] <- value
    return(v)
}
