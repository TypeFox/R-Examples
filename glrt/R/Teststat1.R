Teststat1 <-
function(trt, k=2, est, cens)
{
alpha = t(est$alpha)
m = dim(alpha)[2]
rho = 1-alpha
rho[cens!=3,] = 0
index = est$ppairs
f = est$pf
denom = alpha %*% f

dj = 1:m
for(j in 1:m)
dj[j] = sum(alpha[which(cens!=3),j]*f[j]/denom[which(cens!=3)]) 

djl = matrix(0, nrow=k, ncol=m)   
for(i in 1:k) 
for(j in 1:m) 
djl[i,j] = sum(alpha[which(cens!=3 & trt==i-1),j]*f[j]/denom[which(cens!=3 & trt==i-1)]) 

nj = c(sum(dj), sum(dj)-cumsum(dj)[-m]) + apply(rho, 2, sum) 
njl = matrix(0, nrow=k, ncol=m)
for(i in 1:k)
njl[i,] = c(sum(djl[i,]), sum(djl[i,])-cumsum(djl[i,])[-m]) + apply(rho[which(trt==i-1),], 2, sum) 

u = 1:k
for(i in 1:k)
u[i] = sum(ifelse(nj > 0, djl[i,] - njl[i,] * dj / nj, 0)) 

u
}

