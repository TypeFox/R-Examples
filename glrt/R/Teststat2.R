Teststat2 <-
function(trt, k=2, cens, est, rho=0, gamma=0, c0=1)
{
F = est$sigma
FF = c(0, F)  
index = est$ppairs
tiny = .Machine$double.eps
u = rep(0, k)  
for(i in 1:length(cens))
{
if(cens[i] == 1)  
{
Fu = F[index[i, 2]]
if(Fu > 0 + tiny)  u[trt[i]+1] = u[trt[i]+1] + (Linkfunc(Fu, rho, gamma) - c0)/Fu 
}
else if(cens[i] == 2)  
{
Fu = FF[index[i, 1]]
Fv = F[index[i, 2]]
if(Fv - Fu > 0 + tiny)  u[trt[i]+1] = u[trt[i]+1] + (Linkfunc(Fv, rho, gamma) - Linkfunc(Fu, rho, gamma)) / (Fv - Fu)

}
else  
{
Fv = FF[index[i, 1]]
if(Fv < 1.0 - tiny)  u[trt[i]+1] = u[trt[i]+1] + (c0 - Linkfunc(Fv, rho, gamma)) / (1 - Fv)
}
}
u
}

