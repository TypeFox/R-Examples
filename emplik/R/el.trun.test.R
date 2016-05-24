el.trun.test <- function(y,x,fun=function(t){t},mu,maxit=20,error=1e-9) {
x <- as.vector(x)
y <- as.vector(y)
temp <- Wdataclean2(x,d=rep(1, length(x))) 
x <- temp$value
wt0 <- temp$weight

indi <- function(u,v){ as.numeric(u > v) }
uij <- outer(x,y,FUN="indi")
m <- length(y)
w0 <- rep(1/length(x), length(x))
xmu <- fun(x) - mu
for(i in 1:maxit) {
     pvec0 <- as.vector( w0 %*% uij )
     nvec <- wt0 + as.vector(rowsum( t(w0*(1-uij))/pvec0, group=rep(1, m)))
     w0 <- nvec/sum(nvec)
}
w <- w0
for(i in 1:maxit) {
       pvec <- as.vector( w %*% uij )
       nvec <- wt0 + as.vector(rowsum( t(w*(1-uij))/pvec, group=rep(1, m)))
       w <- el.test.wt(x=xmu, wt=nvec, mu=0)$prob
}
pvec <- as.vector( w %*% uij )
pvec0 <- as.vector( w0 %*% uij )
ELR <- sum(wt0*log(w0)) - sum(log(pvec0)) - sum(wt0*log(w)) + sum(log(pvec))
return(list(NPMLE=w0, NPMLEmu=w, "-2LLR"=2*ELR) )
}
