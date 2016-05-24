if(!require("GNE"))stop("this test requires package GNE.")



#-------------------------------------------------------------------------------
# (1) Example 5 of von Facchinei et al. (2007)
#-------------------------------------------------------------------------------

dimx <- c(1, 1)
#O_i(x)
obj <- function(x, i)
{
	if(i == 1)
		res <- (x[1]-1)^2
	if(i == 2)
		res <- 2*(x[2]-1/2)^2
	res	
}
#Gr_x_j O_i(x)
grobj <- function(x, i, j)
{
	if(i == 1)
		res <- c(2*(x[1]-1), 0)
	if(i == 2)
		res <- c(0, 2*(x[2]-1/2))
	res[j]	
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k)
	2 * (i == j && j == k)

dimlam <- c(1, 1)
#constraint function g_i(x)
h <- function(x)
	sum(x[1:2]) - 1
	
#Gr_x_j h(x)
jach <- function(x)
	matrix(c(1, 1), 1)


x0 <- rexp(sum(dimx))
y <- rexp(sum(dimx))
xy1 <- c(y[1], x0[2])
xy2 <- c(x0[1], y[2])
alpha <- 0.1

resgap <- gapNIR(x0, y, dimx, obj=obj, echo=TRUE)
gapcheck <- obj(x0, 1) - obj(c(y[1], x0[2]), 1) - alpha/2*(x0[1]-y[1])^2 + obj(x0, 2) - obj(c(x0[1], y[2]), 2) - alpha/2*(x0[2]-y[2])^2

tol <- .Machine$double.eps^(1/2)
if(abs(resgap - gapcheck) > tol)
	stop("wrong result 1")

resgrgap <- gradxgapNIR(x0, y, dimx, grobj=grobj) 

grcheck <- c(grobj(x0, 1, 1) - grobj(xy1, 1, 1) + grobj(x0, 2, 1) - grobj(xy2, 2, 1),
grobj(x0, 1, 2) - grobj(xy1, 1, 2) + grobj(x0, 2, 2) - grobj(xy2, 2, 2))
grcheck <- grcheck + c(grobj(xy1, 1, 1),  grobj(xy2, 2, 2))- alpha*(x0-y)


if(sum(abs(resgrgap - grcheck)) > tol)
	stop("wrong result 2")
	

resgrgap <- gradygapNIR(x0, y, dimx, grobj=grobj)

grcheck <- c( - grobj(xy1, 1, 1)  - grobj(xy2, 2, 1),
 - grobj(xy1, 1, 2) - grobj(xy2, 2, 2)) + alpha*(x0-y)


if(sum(abs(resgrgap - grcheck)) > tol)
	stop("wrong result 3")
	
	 
	
	
fpNIR(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, echo=FALSE, control=list(eps=1e-10))	
fpNIR(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, echo=FALSE, yinit=runif(2, max=1/2))



fpNIR(x0, dimx, obj=obj, grobj=grobj, echo=FALSE)		
fpNIR(x0, dimx, obj=obj, joint=h, echo=FALSE)		
fpNIR(x0, dimx, obj=obj, echo=FALSE)		

res1 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="pure", merit="FP", control.outer=list(maxiter=10, trace=1))
res1$inner.counts
res1$outer.counts


res2 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="pure", merit="NI", control.outer=list(maxit=10, echo=2))
res2$inner.counts
res2$outer.counts

if(FALSE)
{
res3 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="UR", merit="FP", control.outer=list(maxit=10), stepfunc=decrstep, argstep=5)

res3 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="UR", merit="FP", control.outer=list(maxit=10, echo=3), stepfunc=decrstep5)
res3$inner.counts
res3$outer.counts

res4 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="UR", merit="NI", control.outer=list(maxit=10, echo=3), stepfunc=decrstep5)
res4$inner.counts
res4$outer.counts

res5 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="RRE", merit="NI", control.outer=list(maxit=10, echo=3))
res5$inner.counts
res5$outer.counts

res5 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="RRE", merit="FP", control.outer=list(maxiter=10, trace=1), order=1)
res5$inner.counts
res5$outer.counts

res6 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="MPE", merit="NI", control.outer=list(maxit=10, echo=3))
res6$inner.counts
res6$outer.counts

res6 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="MPE", merit="FP", control.outer=list(maxiter=10, trace=1), order=1)
res6$inner.counts
res6$outer.counts

res7 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="SqMPE", merit="NI", control.outer=list(maxit=10, echo=3))
res7$inner.counts
res7$outer.counts

res7 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="SqMPE", merit="FP", control.outer=list(maxiter=10, trace=1), order=1)
res7$inner.counts
res7$outer.counts

res8 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="SqRRE", merit="NI", control.outer=list(maxit=10, echo=3))
res8$inner.counts
res8$outer.counts

res8 <- GNE.fpeq(x0, dimx, obj=obj, grobj=grobj, joint=h, jacjoint=jach, method="SqRRE", merit="FP", control.outer=list(maxiter=10, trace=1))
res8$inner.counts
res8$outer.counts
}


#-------------------------------------------------------------------------------
# (5) another Example
#-------------------------------------------------------------------------------

# associated objective functions

dimx <- c(2, 2, 3)
#O_i(x)
fullob <- function(x, i)
{
	x <- x[1:7]	
	if(i == 1)
		res <- sum((x - 1:7)^3)
	if(i == 2)
		res <- sum((x - 1:7)^(1:7))
	if(i == 3)
		res <- x[1] + x[3] + (x[5]^2 + x[6]^2 + x[7]^2 - 5)^2
	
	res
}
#Gr_x_j O_i(x)
grfullob <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
		grad <- 3*(x - 1:7)^2
	if(i == 2)
		grad <- 1:7*(x - 1:7)^(0:6)
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2 - 5	
		grad <- c(1, 0, 1, 0, 4*x[5]*s, 4*x[6]*s, 4*x[7]*s)
	}
	grad[j]	
}


#Gr_x_k Gr_x_j O_i(x)
hefullob <- function(x, i, j, k)
{
	x <- x[1:7]
	if(i == 1)
		he <- diag( 6*(x - 1:7) )
	if(i == 2)
		he <- diag( c(0, 2, 6, 12, 20, 30, 42)*(x - 1:7)^c(0, 0:5) )
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2	
		
		he <- rbind(rep(0, 7), rep(0, 7), rep(0, 7), rep(0, 7),
					c(0, 0, 0, 0, 4*s+8*x[5]^2, 8*x[5]*x[6], 8*x[5]*x[7]),
					c(0, 0, 0, 0, 8*x[5]*x[6], 4*s+8*x[6]^2, 8*x[6]*x[7]),
					c(0, 0, 0, 0,  8*x[5]*x[7], 8*x[6]*x[7], 4*s+8*x[7]^2) )
	}
	he[j,k]	
}

x0 <- rexp(sum(dimx))
y <- rexp(sum(dimx))
xy1 <- c(y[1:2], x0[3:7])
xy2 <- c(x0[1:2], y[3:4], x0[5:7])
xy3 <- c(x0[1:4], y[5:7])
alpha <- 0.1

resgap <- gapNIR(x0, y, dimx, obj=fullob)

gapcheck <- fullob(x0, 1)-fullob(xy1, 1)-alpha/2*sum((x0[1:2]-y[1:2])^2)
gapcheck <- gapcheck + fullob(x0, 2)-fullob(xy2, 2)-alpha/2*sum((x0[3:4]-y[3:4])^2)
gapcheck <- gapcheck + fullob(x0, 3)-fullob(xy3, 3)-alpha/2*sum((x0[5:7]-y[5:7])^2)

if(sum(abs(resgap - gapcheck)) > tol)
	stop("wrong result 4")


resgrxgap <- gradxgapNIR(x0, y, dimx, grobj=grfullob)

grxcheck <- sapply(1:7, function(j) grfullob(x0, 1, j)) - sapply(1:7, function(j) grfullob(xy1, 1, j)) 
grxcheck <- grxcheck + sapply(1:7, function(j) grfullob(x0, 2, j)) - sapply(1:7, function(j) grfullob(xy2, 2, j)) 
grxcheck <- grxcheck + sapply(1:7, function(j) grfullob(x0, 3, j)) - sapply(1:7, function(j) grfullob(xy3, 3, j)) 
grxcheck <- grxcheck + c(sapply(1:2, function(j) grfullob(xy1, 1, j)), sapply(3:4, function(j) grfullob(xy2, 2, j)), sapply(5:7, function(j) grfullob(xy3, 3, j)))
grxcheck <- grxcheck - alpha*(x0-y) 

if(sum(abs(resgrxgap - grxcheck)) > tol)
	stop("wrong result 5")


resgrygap <- gradygapNIR(x0, y, dimx, grobj=grfullob)

grycheck <- - c(sapply(1:2, function(j) grfullob(xy1, 1, j)), sapply(3:4, function(j) grfullob(xy2, 2, j)), sapply(5:7, function(j) grfullob(xy3, 3, j)))
grycheck <- grycheck + alpha*(x0-y) 

if(sum(abs(resgrygap - grycheck)) > tol)
	stop("wrong result 6")


