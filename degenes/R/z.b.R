z.b <-
function (a, b, c, d, m1, m2) 
{
    cat("Calculation of Z", fill = TRUE)
    Z <- (apply(a, 1, mean) - apply(b, 1, mean))/sqrt((apply(c, 
        1, var)/m1) + (apply(d, 1, var)/m1))
    return(Z)
}

