loss <-
function(par,x)
{
p <- par[1]
q <- par[2]
M <- Edlaplace2(p, q)
(mean(x)-M$E1)^2+(mean(x^2)-M$E2)^2
}
