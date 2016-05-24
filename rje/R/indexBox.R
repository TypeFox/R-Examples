indexBox <-
function (upp, lwr, dim) 
{
    if (any(upp < 0) || any(lwr > 0)) 
        stop("Incorrect bounds")
    if (any(dim <= 0)) 
        stop("Incorrect dimensions")
    ld = length(dim)
    if (length(upp) != ld) 
        upp = rep(upp, length.out = ld)
    if (length(lwr) != ld) 
        lwr = rep(lwr, length.out = ld)
    out = .C("indexBox", as.integer(upp), as.integer(lwr), as.integer(dim), 
        ld, integer(prod(upp - lwr + 1)), package = "rje")
    return(out[[5]])
}
