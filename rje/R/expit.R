expit <-
function (x) 
{
    out = exp(x)/(1 + exp(x))
    out[x > 100] = 1
    out
}
