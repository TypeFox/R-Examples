ci.c.ancova<-function(Psi=NULL, adj.means=NULL, s.ancova=NULL, c.weights, n, cov.means, SSwithin.x, conf.level=.95, ...)
{
if (length(cov.means)!=length(c.weights) ) stop("The input 'cov.means' and 'c.weights' imply different number of groups")

if(is.null(Psi) & is.null(adj.means) ) stop("Input either 'Psi' or 'adj.means'")
if(!is.null(Psi) & !is.null(adj.means) ) stop("Do not input both 'Psi' and 'adj.means' at the same time")

if(is.null(Psi)) Psi<- sum(adj.means*c.weights)
J<- length(c.weights)
if(length(n)==1) n<-rep(n, J)
if(length(n)>1 & length(n)!=length(c.weights)) stop("The input 'n' and 'c.weights' imply different number of groups ")
########################################################################
f.x.numerater<- ( sum(c.weights*cov.means) )^2
f.x.denominator<- SSwithin.x
sample.size.weighted<- sum(c.weights^2 / n)

se.Psi<- s.ancova*sqrt(sample.size.weighted + f.x.numerater/f.x.denominator)

alpha<- 1-conf.level
nu<- sum(n)-J-1    
t.value<- qt(1-alpha/2, df=nu) 

list(lower.limit=Psi - t.value*se.Psi, Psi=Psi, upper.limit=Psi + t.value*se.Psi)
}
