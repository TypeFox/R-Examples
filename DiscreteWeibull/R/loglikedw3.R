loglikedw3<-function(par,x)
{
sum(-log(ddweibull3(x,par[1],par[2])))
}