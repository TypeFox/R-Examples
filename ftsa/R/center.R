center = function (ab, k) 
{
    h <- diff(ab)/k
    list(seq(ab[1] + h/2, by = h, length = k))
}
