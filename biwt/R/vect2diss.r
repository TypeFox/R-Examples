vect2diss <- 
function (v) 
# This function is very similar to 
# dissmatrix in the hopach package
# except that it fills in by the lower
# triangle of the matrix instead of the
# upper triangle

{
    if (!is.vector(v)) 
        stop("arg to dissmatrix() must be a vector")
    p <- (1 + sqrt(1 + 8 * length(v)))/2
    M <- matrix(0, nrow = p, ncol = p)
    count <- 1
    for (i in 1:(p - 1)) {
        M[(i+1), 1:i] <- v[count:(count + i - 1)]
        count <- count + i
    }
    return(M + t(M))
}