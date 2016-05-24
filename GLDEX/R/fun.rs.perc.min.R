"fun.rs.perc.min" <-
function(coef, u = 0.1, data)
{
L3 <- coef[1]
L4 <- coef[2]
pp <- fun.lambda.percentile(data)
rho3 <- ((1 - u)^L4 - u^L3 + 0.5^L3 - 0.5^L4)/((1 - u)^L3 - u^L4 + 0.5^L4 - 0.5^L3)
rho4 <- (0.75^L3 - 0.25^L4 + 0.75^L4 - 0.25^L3)/((1 - u)^L3 - u^L4 + (1 - u)^L4 - u^L3)
sqrt((rho3 - pp$p3)^2 + (rho4 - pp$p4)^2)
}

