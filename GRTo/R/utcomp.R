utcomp <-
function (b1, n1, b2, n2) 
{
    da = -2 * (n1 + n2) * log(n1 + n2) + 2 * n1 * log((n2 * b1)/b2 + 
        n1) + 2 * n2 * log((n1 * b2)/b1 + n2) - 2
    p_ut = exp(-da/2 - 2)
    p_ut
}
