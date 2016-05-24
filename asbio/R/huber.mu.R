huber.mu<-function(x,c=1.28,iter=20,conv=.1e-06){
mu.hat <- huber.NR(x, c, iter)
if(c < 0) stop("c must be >= 0")
if(c == 0 | is.nan(mu.hat[2])) mu.est <- median(x)
if(c > 0 & !is.nan(mu.hat[2]))
for(i in 1:iter){
if(abs(mu.hat[i] - mu.hat[i+1]) <=conv) mu.est=mu.hat[i+1]
}
mu.est
}
