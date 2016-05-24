Teststat3 <-
function(trt, k, cens, counts, est, rho=0, gamma=0, c0=1)
{

N = apply(counts, 2, sum)  
F = est$sigma
FF = c(0, F)  
index = est$ppairs
u = rep(0, k)  
tiny = .Machine$double.eps
for(i in 1:length(cens))
{


if(cens[i] == 1)  
{
Fu = F[index[i, 2]]
if(Fu > 0 + tiny)
u[trt[i]+1] = u[trt[i]+1] + (Linkfunc(Fu, rho, gamma) - c0)/Fu * N[2] / counts[trt[i]+1,2]
}
else if(cens[i] == 2)  
{
Fu = FF[index[i, 1]]
Fv = F[index[i, 2]]
if(Fv - Fu > 0 + tiny)
u[trt[i]+1] = u[trt[i]+1] + (Linkfunc(Fv, rho, gamma) - Linkfunc(Fu, rho, gamma)) / (Fv - Fu) * N[2] / counts[trt[i]+1,2]
}
else if(cens[i] == 3)  
{
Fv = FF[index[i, 1]]
if(Fv < 1.0 - tiny)
u[trt[i]+1] = u[trt[i]+1] + (c0 - Linkfunc(Fv, rho, gamma)) / (1 - Fv) * N[2] / counts[trt[i]+1,2]
}

else
{
Ft = FF[index[i,2]+1]

Ftminus = FF[index[i,2]]

if(abs(Ft - Ftminus) > tiny)  # Ft != Ftminus

u[trt[i]+1] = u[trt[i]+1] + ( (Linkfunc(Ft, rho, gamma) - Linkfunc(Ftminus, rho, gamma)) / (Ft - Ftminus) ) * N[1] / counts[trt[i]+1,1]

}
}
u
}

