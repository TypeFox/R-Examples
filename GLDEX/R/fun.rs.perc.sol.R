"fun.rs.perc.sol" <-
function(coef, u = 0.1, data)
{
L3 <- coef[1]
L4 <- coef[2]
pp <- fun.lambda.percentile(data)
L2 <- ((1 - u)^L3 - u^L4 + (1 - u)^L4 - u^L3)/pp$p2
L1 <- pp$p1 - (0.5^L3 - 0.5^L4)/L2
return(list("L1"=L1, "L2"=L2, "L3"=L3, "L4"=L4))
}

