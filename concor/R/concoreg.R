"concoreg" <-
function(x,y,py,r) { 
n<-dim(x)[1]
q<-dim(y)[2]
if (sum(py) != q ) stop("py IS NOT SUITABLE")

# rk here is a maximal value given for the rank of x. 
# you may modify the tolerance 1e-8
s<-svd(x);rk<-sum(s$d > max(dim(x))*s$d[1]*1e-8)
if (r > min(c(min(py),rk,n))) stop("r IS TOO HIGH")

P=matrix(s$u[,1:rk]*sqrt(n),ncol=rk)
s<-concor(P,y,py,r)
list(cx=P%*%s$u,v=s$v,V=s$V,varexp=s$cov2)
}

