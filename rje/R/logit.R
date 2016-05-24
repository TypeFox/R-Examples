logit <-
function (x) 
{
    no = (x < 0) | (x > 1)
    out = numeric(length(x))
    out[no] = NaN
    out[!no] = log(x[!no]/(1 - x[!no]))
    dim(out) = dim(x)
    out
}
