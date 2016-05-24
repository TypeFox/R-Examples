ScoreStat <-
function(trt, k, est)
{
alpha = t(est$alpha)
n = nrow(alpha)
m = ncol(alpha) 
u = rep(0, k)  
Ghat = c(1, 1-est$sigma) 
alpha = cbind(rep(0, n), alpha)  
ghat = -diff(Ghat)   
phat = rep(0, m)  
tiny = .Machine$double.eps   
for(j in 1:m) 
{ 
if(Ghat[j] > tiny)  # Ghat[j] != 0
phat[j] = Ghat[j+1] / Ghat[j]
}
denom = alpha[,-1] %*% ghat
for(i in 1:n)
{
num = alpha[i,-1] * ghat
for(j in 1:m)
{
if(abs(phat[j] - 1) > 0+tiny && abs(phat[j]) > 0+tiny)
u[trt[i]+1] = u[trt[i]+1] + log(phat[j]) * sum(num[j:m])/ denom[i] - log(phat[j]) * alpha[i,j+1] * ghat[j] / (1-phat[j]) / denom[i]
}
}
u
}

