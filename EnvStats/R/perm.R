perm <-
function (J, n) 
{
    if (!is.vector(J, mode = "numeric") || is.factor(J) || length(J) != 
        1 || J != trunc(J)) 
        stop("'J' must be an integer")
    if (!is.vector(n, mode = "numeric") || is.factor(n) || length(n) != 
        1 || n != trunc(n)) 
        stop("'n' must be an integer")
    if (J < 1 || J > factorial(n)) 
        stop("'J' must be an integer between 1 and 'n'")
    if (n == 0) 
        return(integer(0))
    if (n == 1) 
        return(as.integer(1))
    J <- J - 1
    per <- 1
    for (ind in 2:n) {
        nleft <- J%%ind
        J <- J%/%ind
        if (nleft == 0) 
            per <- c(ind, per)
        else {
            if (nleft == ind - 1) 
                per <- c(per, ind)
            else per <- c(per[1:nleft], ind, per[(1 + nleft):(ind - 
                1)])
        }
    }
    per
}
