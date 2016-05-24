solbeta=function(a,b,c,prec=10^(-3))
{
coc=(1-c)/c
detail=alpha=1
while (probet(a,b,c,alpha)<.95) alpha=alpha+detail
while (abs(probet(a,b,c,alpha)-.95)>prec)
{
alpha=max(alpha-detail,detail/10)
detail=detail/10
while (probet(a,b,c,alpha)<.95) alpha=alpha+detail
}
list(alpha=alpha,beta=alpha*coc)
}
