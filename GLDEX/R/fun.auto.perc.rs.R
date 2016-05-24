"fun.auto.perc.rs" <-
function(data, init, u = 0.1, method = "Nelder-Mead")
{
# Find solution for L3 and L4
optim.r <- optim(par = init, gr = fun.rs.perc.gradient, 
fn = fun.rs.perc.min, data = data, method = method)
p.sol <- optim.r$par
value <- optim.r$value
msg <- optim.r$message
# Find all the solutions:
all.sol <- fun.rs.perc.sol(p.sol, data = data)
return(list("all.sol"=all.sol, "value"=value, "msg"=msg))
}

