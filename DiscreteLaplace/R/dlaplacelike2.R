dlaplacelike2 <-
function(par,x)
{
p <- par[1]
q <- par[2]
sum(-log(ddlaplace2(x,p,q)))
}
