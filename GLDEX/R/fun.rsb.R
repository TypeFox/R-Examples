"fun.rsb" <-
function(L3, L4, k)
{
j <- 0:k
result <- sum(((-1)^j * (gamma(k + 1)/(gamma(k - j + 1) * gamma(j + 1))) * fun.beta(L3 * (k - j) + 1,L4 * j + 1)))
}

