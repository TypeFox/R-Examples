"fun.fmkl.mm.sol.alt" <-
function(result, data)
{
L3 <- result[,1]
L4 <- result[,2]
aa <- fun.moments(data)
v1 <- (L3 * (L3 + 1))^-1 - (L4 * (L4 + 1))^-1
v2 <- (L3^2 * (2 * L3 + 1))^-1 + (L4^2 * (2 * L4 + 1))^-1 - (2/(L3 * L4)) * fun.beta(L3 + 1, L4 + 1)
L2 <- sqrt((v2 - v1^2)/aa$a2)
L1 <- aa$a1 + (1/L2) * ((L3 + 1)^-1 - (L4 + 1)^-1)
return(list("L1"=L1, "L2"=L2, "L3"=L3, "L4"=L4))
}

