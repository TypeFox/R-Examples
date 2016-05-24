Z2freq <-
function (Z) 
{
    s <- apply(Z, 2, sum)
    return(s/(sum(s)))
}
