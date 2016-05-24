#
# This demo illustrates that the imprecise estimate found from our
# approach is identical to those found from constrained optimization
# technique.  
# 
# updated on 2015.11.23
# 

#
# case 1. imprecise log-gamma prior (2 dims)
#

# characterization of imprecise prior 
ui <- rbind(diag(2), -diag(2))
ci <- c(0,0,-1,-1)
op <- iprior(ui=ui, ci=ci)
op1 <- update(op, y=1)
sop1 <- summary(op1)

# constrained optimization 
ui <- rbind(diag(2), -diag(2))
ci <- c(1,1,-2,-2)
m0 <- stats::constrOptim(theta=c(1.5, 1.5), f=fn.evfn, grad=gr.evfn, ui=ui, ci=ci, method="BFGS", hessian=TRUE)
m1 <- stats::constrOptim(theta=c(1.5, 1.5), f=fn.evfn, grad=gr.evfn, ui=ui, ci=ci, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

# optima found by constrained optimization
m0$par
m1$par

# optima found by our approach
sop1$inf.p1
sop1$sup.p1

# lower and upper bounds found by constrained optimization 
c(lower=m0$value, upper=m1$value)

# lower and upper bounds found by our approach 
c(lower=sop1$inf, upper=sop1$sup)


#
# case 2. imprecise normal prior (3 dim)
#

# characterization of imprecise prior 
ui <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(-1,0,0),c(0,-1,0),c(0,0,-1))
ci <- c(1,1,1,-2,-2,-2)
op <- iprior(ui=ui, ci=ci)

# constrained optimization 
m0 <- stats::constrOptim(theta=c(1.5,1.5,1.5), f=fn.evfn, grad=gr.evfn, ui=ui, ci=ci, method="BFGS", hessian=TRUE)
m1 <- stats::constrOptim(theta=c(1.5,1.5,1.5), f=fn.evfn, grad=gr.evfn, ui=ui, ci=ci, method="BFGS", hessian=TRUE, control=list(fnscale=-1))

# our approach
op1 <- iprior(pmat=op$vtx)
op1 <- update(op1, y=numeric(0))
sop1 <- summary(op1)

# optima found by constrained optimization
m0$par
m1$par

# optima found by our approach
sop1$inf.p1
sop1$sup.p1

# lower and upper bounds found by constrained optimization 
c(lower=m0$value, upper=m1$value)

# lower and upper bounds found by our approach 
c(lower=sop1$inf, upper=sop1$sup)

