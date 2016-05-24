lossdiw <-
function(x, par, eps=0.0001, nmax=1000)
{
E<-Ediweibull(par[1],par[2],eps=eps, nmax=nmax)

EX<-E[[1]]
EX2<-E[[2]]

(mean(x)-EX)^2+(mean(x^2)-EX2)^2
}

