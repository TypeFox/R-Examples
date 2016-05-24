lossdw <-
function(par,x,zero=FALSE,eps=0.0001, nmax=1000)
{(mean(x)-Edweibull(par[1],par[2],zero=zero,eps=eps,nmax=nmax))^2+(mean(x^2)-E2dweibull(par[1],par[2],zero=zero,eps=eps,nmax=nmax))^2}

