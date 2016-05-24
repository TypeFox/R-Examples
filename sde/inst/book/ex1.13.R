# ex1.13.R
set.seed(123)
n <- 5
N <- n*10
t <- seq(0, 2*pi, length=N+1)
f <- sin(t)+rnorm(N+1)
plot(t, f, type="l", axes=F, 
   ylab=expression(f(t)==sin(t)+epsilon[t]))
idx <- seq(1, N+1, length=n+1)
axis(1,t[idx], c(sprintf("%3.2f", t[idx[1:n]]),
   expression(2*pi)))
axis(2)
box()
for(i in 1:n){
   lines( c(t[idx[i]],t[idx[i+1]]), c(f[idx[i]],f[idx[i]]) )
   points( t[idx[i]], f[idx[i]] ,pch=19)
}
text(t[idx[3]]+.5,f[idx[2]]+.5,expression(I(f^{(5)})))
