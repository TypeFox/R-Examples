gl.unequal <-
function (n = "number of levels", k = "numeric vector of sample sizes") 
{
    out <- numeric()
    for (i in 1:n) {
        out <- append(out, rep(i, length.out = k[i]))
    }
    return(factor(out))
}
