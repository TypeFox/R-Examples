"fun.rs.nr" <-
function(x, l1, l2, l3, l4, init, error = 1e-010)
{
u <- init
while(abs(l1 + (u^l3 - (1 - u)^l4)/l2 - x) > error) {
u <- u - ( - x + l1 + (u^l3 - (1 - u)^l4)/l2)/((l3 * u^(l3 - 1) + l4 * (1 - u)^(l4 - 1))/l2)
u <- u
}
return(u)
}

