"ifa.init.random" <-
function(y,L)
{
numvar<-ncol(y)
H<-matrix(runif(L*numvar,-3,3),numvar,L)
psi<-diag(diag(var(y)))*0.1

output<-list(psi=psi,H=H)
}
