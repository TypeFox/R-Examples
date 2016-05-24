cpo <-
function (obj) 
{
    fx <- obj$fx
    cpo <- 1/apply(1/fx, 1, mean)
    return(cpo)
}
