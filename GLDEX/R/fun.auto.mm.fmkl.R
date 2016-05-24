"fun.auto.mm.fmkl" <-
function(data, init, method = "Nelder-Mead")
{
# Find solution for L3 and L4
optim.r <- optim(init, fun.fmkl.mm.min, data = data, method = method)
p.sol <- optim.r$par
value <- optim.r$value
msg <- optim.r$message
# Find all the solutions:
all.sol <- fun.fmkl.mm.sol(p.sol, data)
return(list("all.sol"=all.sol, "value"=value, "msg"=msg))
}

