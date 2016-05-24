VarScoreStat <-
function(trt, k, est)
{
alpha = t(est$alpha)
n = nrow(alpha)
m = ncol(alpha)
Ghat = c(1, 1-est$sigma) 
alpha = cbind(rep(0, n), alpha)  
ghat = -diff(Ghat)   

b = rep(0, m+1)
h = rep(0, m+1)
tiny = .Machine$double.eps
for(j in 1:(m+1))
{
if(Ghat[j] > 0 + tiny)
{
h[j] = Ghat[j] * log(Ghat[j])
b[j] = Ghat[j] * log(Ghat[j]) + Ghat[j] * log(Ghat[j])^2
} 
}

denom = alpha[,-1] %*% ghat  
alphadiff = alpha[,3:(m+1)] - alpha[,2:m]

mu = alphadiff/matrix(rep(denom,m-1), ncol=m-1)

A11 = matrix(0, nrow=m-1, ncol=m-1)
for(j in 1:(m-1))
{
for(r in 1:(m-1))
{
if(j != r)
A11[j,r] = A11[j,r] + sum(mu[,j] * mu[,r] * h[j+1] * h[r+1])
else
A11[j,r] = A11[j,r] + sum(-mu[,j] * b[j+1] + (mu[,j] * h[j+1])^2)
}
}

A12 = matrix(0, nrow=k, ncol=m-1)
sum12 = alpha[,2:(m+1)] %*% (-diff(h))
for(i in 1:n)
{
for(r in 1:k)
{
for(j in 1:(m-1))
{
if(r == trt[i] + 1) 
{
A12[r, j] = A12[r, j] - mu[i, j] * b[j+1] + mu[i,j] * h[j+1] * sum12[i] / denom[i] 
}
}
}
}


A22 = matrix(0, nrow=k, ncol=k)
sum221 = alpha[,2:(m+1)] %*% (-diff(b))
sum222 = alpha[,2:(m+1)] %*% (-diff(h))
for(i in 1:n)
{
for(j in 1:k)  
for(r in 1:k)
if(j == trt[i] + 1 && r == trt[i] + 1)
A22[j, r] = A22[j, r] - sum221[i] / denom[i] + (sum222[i] / denom[i])^2  
}
V = A22 - A12 %*% solve(A11) %*% t(A12)
V
}

