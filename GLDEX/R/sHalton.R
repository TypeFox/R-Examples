"sHalton" <-
function (n.max, n.min = 1, base = 2, leap = 1) 
{
    stopifnot((leap <- as.integer(leap)) >= 1)
    nd <- as.integer(1 + log(n.max, base))
    dB <- digitsBase(if (leap == 1) 
        n.min:n.max
    else seq(n.min, n.max, by = leap), base = base, ndigits = nd)
    colSums(dB/base^(nd:1))
}

