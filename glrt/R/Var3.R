Var3 <-
function(trt, k, cens, counts, est, rho=0, gamma=0, c0=1)
{
N = apply(counts, 2, sum)
n = length(trt)
nn = 1:k;
for(i in 1:k) nn[i] = length(trt[trt==i-1])
F = est$sigma
FF = c(0, F)  
index = est$ppairs
Q = 0
P = 0
tiny = .Machine$double.eps
for(i in 1:n)
{
if(cens[i] == 1)  
{
Fu = F[index[i, 2]]
if(Fu > tiny)  Q = Q + ( (Linkfunc(Fu, rho, gamma) - c0)/Fu )**2 / N[2]
}
else if(cens[i] == 2)  
{
Fu = FF[index[i, 1]]
Fv = F[index[i, 2]]
if(Fv - Fu > tiny)
Q = Q + ( (Linkfunc(Fv, rho, gamma) - Linkfunc(Fu, rho, gamma)) / (Fv - Fu) )**2 / N[2]
}
else if(cens[i] == 3)  
{
Fv = FF[index[i, 1]]
if(Fv < 1.0 - tiny)
Q = Q + ( (c0 - Linkfunc(Fv, rho, gamma)) / (1 - Fv) )**2 / N[2]
}
else  
{
Fv = F[index[i,2]]
Fvminus = FF[index[i,1]]

if(abs(Fv - Fvminus) <= tiny)  # Fv == Fvminus
P = P + ( LinkfuncDir(Fv, rho, gamma) )^2 / N[1]
else
P = P + ( (Linkfunc(Fv, rho, gamma) - Linkfunc(Fvminus, rho, gamma)) / (Fv - Fvminus) )^2 / N[1]
}
}
v = matrix(0, nrow=k, ncol=k)  
for(i in 1:k)
{
for(j in 1:k)
{
if(i == j)
{
if(counts[j,1] != 0)
v[i,j] = P * N[1] * (N[1]/counts[j, 1] - 1) / n
v[i,j] = v[i,j] + Q * N[2] * (N[2]/counts[j,2] - 1) / n
}
else
v[i,j] = -P * N[1] / n - Q * N[2] / n
}
}
v
}

