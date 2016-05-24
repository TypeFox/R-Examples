combns <-
function (x) 
{
    ans <- t(combn(seq(x), 2))
    ans <- cbind(ans, seq.int(nrow(ans)))
    ans
}

