Var1 <-
function(A, est, trt, cens, k=2, M=50)
{
f = est$pf
index = est$ppairs
intmap = est$intmap
A2 = A  

u = matrix(0, nrow=M, ncol=k)
var = matrix(0, nrow=k, ncol=k)
cens2 = cens
for(r in 1:M)
{
for(i in 1:length(cens))  
{
if(cens[i]!=3)
{
F = cumsum( f[index[i,1]:index[i,2]]/sum(f[index[i,1]:index[i,2]]) ) 
t = intmap[2, index[i,1]:index[i,2]] 
A2[i, 2] = t[length(F[F < runif(1)]) + 1]
 A2[i, 1] = A2[i,2] 
cens2[i] = 1
}
else  cens2[i] = 0
}

est2 = survdiff(Surv(A2[,1], cens2) ~ trt) 
u[r,] = est2$obs - est2$exp
var = var + est2$var
}
v1 = var / M
Umean = apply(u, 2, mean)
v2 = matrix(0, nrow=k, ncol=k)
for(r in 1:M)
v2 = v2 + outer(u[r,]-Umean, u[r,]-Umean)
v2 = (1 + 1/M) * v2 / (M-1)
v = v1 + v2
v
}

