Var2 <-
function(trt, k=2, cens, est, rho=0, gamma=0, c0=1)
{
n = length(trt)
nn = 1:k;
for(i in 1:k) nn[i] = length(trt[trt==i-1])

F = est$sigma
FF = c(0, F)  
index = est$ppairs

Q = 0
tiny = .Machine$double.eps
for(i in 1:length(cens))
{
if(cens[i] == 1)  
{
Fu = F[index[i, 2]]
if(Fu > 0 + tiny)  Q = Q + ( (Linkfunc(Fu, rho, gamma) - c0)/Fu )**2 
}
else if(cens[i] == 2)  
{
Fu = FF[index[i, 1]]
Fv = F[index[i, 2]]
if(Fv - Fu > 0 + tiny)  Q = Q + ( (Linkfunc(Fv, rho, gamma) - Linkfunc(Fu, rho, gamma)) / (Fv - Fu) )**2

}
else  
{
Fv = FF[index[i, 1]]
if(Fv < 1.0 - tiny)  Q = Q + ( (c0 - Linkfunc(Fv, rho, gamma)) / (1 - Fv) )**2
}
}

v = matrix(0, nrow=k, ncol=k)  
for(i in 1:k)
{
for(j in 1:k)
{
if(i == j)  v[i,j] = Q * nn[i] * (n - nn[i]) / n**2
else  v[i,j] = -Q * nn[i] * nn[j] / n**2
}
}
v
}

