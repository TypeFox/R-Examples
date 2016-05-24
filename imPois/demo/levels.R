# 
# This demo shows a set of level surface from -3 to 3 running by 1 
# on the expectation of canonical parameter. 
# 
# updated on 2015.11.21
# 

# Grid on the hyperparameter space 
e2 <- seq(from=0.05, to=4, by=0.25)
e0 <- seq(from=0, to=2.5, by=0.5)
et <- seq(from=-3, to=3, by=1)
h0 <- expand.grid(e2=e2, e0=e0, et=et)
fn <- function(x){
	f0 <- function(xm, ...) evfn(y=numeric(0), pars=c(x[1], xm, x[2]))$value - x[3]
	val <- tryCatch(uniroot(f0, lower=-10, upper=10, extendInt="yes", tol=1e-5)$root, error=function(e) return(NaN))
	return(val)	
}
h0$e1 <- apply(h0, 1, fn)

# Set of nonlinear level surfaces
require(rgl)
m <- matrix(h0$e1, ncol=length(e0), nrow=length(e2))
rgl.open()
rgl.clear()
rgl.bg(color="white")
rgl.points(x=h0$e0, y=h0$e2, z=h0$e1, color="grey")
decorate3d()
aspect3d(1,1.5,2)

cols <- rainbow(7, s=0.7)
for(i in -3:3){
	h1 <- subset(h0, subset=(et==i))
	m1 <- matrix(h1$e1, ncol=length(e0), nrow=length(e2))
	surface3d(x=e0, y=e2, z=t(m1), color=cols[i+4])
}

