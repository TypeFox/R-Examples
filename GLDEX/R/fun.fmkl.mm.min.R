"fun.fmkl.mm.min" <-
function(coef, data)
{
L3 <- coef[1]
L4 <- coef[2]
aa <- fun.moments(data)
v1 <- fun.fmklb(L3, L4, 1)
v2 <- fun.fmklb(L3, L4, 2)
v3 <- fun.fmklb(L3, L4, 3)
v4 <- fun.fmklb(L3, L4, 4)
g3 <- (v3 - 3 * v2 * v1 + 2 * v1^3) * (v2 - v1^2)^(-3/2)
g4 <- (v4 - 4 * v1 * v3 + 6 * v1^2 * v2 - 3 * v1^4) * (v2 - v1^2)^(-2)
abs(g3 - aa$a3) + abs(g4 - aa$a4)
}

