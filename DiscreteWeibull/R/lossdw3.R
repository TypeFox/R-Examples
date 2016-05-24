lossdw3<-function(par,x,eps=1e-04)
{
(mean(x)-Edweibull3(par[1],par[2],eps=eps))^2+(mean(x^2)-E2dweibull3(par[1],par[2],eps=eps))^2
}