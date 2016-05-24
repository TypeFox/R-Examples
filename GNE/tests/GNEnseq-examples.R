if(!require("GNE"))stop("this test requires package GNE.")

itermax <- 10

#-------------------------------------------------------------------------------
# (1) Example 5 of von Facchinei et al. (2007)
#-------------------------------------------------------------------------------

dimx <- c(1, 1)
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
g <- function(x, i)
	sum(x[1:2]) - 1
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	1
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0



#true value is (3/4, 1/4, 1/2, 1/2)

z0 <- rep(0, sum(dimx)+sum(dimlam))

funSSR(z0, dimx, dimlam, grobj=grobj, constr=g, grconstr=grg, compl=phiFB, echo=FALSE)

	
jacSSR(z0, dimx, dimlam, heobj=heobj, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, NULL, heobj=heobj, NULL, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1, maxit=itermax))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, NULL, heobj=heobj, NULL, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1, maxit=itermax))



#-------------------------------------------------------------------------------
# (2) Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------


#constants
myarg <- list(d= 20, lambda= 4, rho= 1)

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	res <- -arg$rho * x[i]
	if(i == j)
		res <- res + arg$d - arg$lambda - arg$rho*(x[1]+x[2])
	-res
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
	arg$rho * (i == j) + arg$rho * (j == k)	


dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	-x[i]
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	-1*(i == j)
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0

#true value is (16/3, 16/3, 0, 0) 

z0 <- rep(0, sum(dimx)+sum(dimlam))

funSSR(z0, dimx, dimlam, grobj=grobj, myarg, constr=g, grconstr=grg, compl=phiFB, echo=FALSE)

jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1, maxit=itermax))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1, maxit=itermax))


#-------------------------------------------------------------------------------
# (3) River basin pollution game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------

myarg <- list(
C = cbind(c(.1, .12, .15), c(.01, .05, .01)),
U = cbind(c(6.5, 5, 5.5), c(4.583, 6.25, 3.75)),
K = c(100, 100),
E = c(.5, .25, .75),
D = c(3, .01)
)



dimx <- c(1, 1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	dij <- 1*(i == j)
	res <- -(-arg$D[2] - arg$C[i, 2]*dij) * x[i] 
	res - (arg$D[1] - arg$D[2]*sum(x[1:3]) - arg$C[i, 1] - arg$C[i, 2]*x[i]) * dij
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
	dij <- 1*(i == j)
	dik <- 1*(i == k)
	
	arg$D[2] * dik + arg$D[2] * dij + 2 * arg$C[i, 2] * dij * dik
}

dimlam <- c(2, 2, 2)
#g_i(x)
g <- function(x, i, arg)
	c(sum(arg$U[, 1] * arg$E * x[1:3]) - arg$K[1],
	  sum(arg$U[, 2] * arg$E * x[1:3]) - arg$K[2])
#Gr_x_j g_i(x)
grg <- function(x, i, j, arg)
	c(arg$U[j, 1] * arg$E[j], arg$U[j, 2] * arg$E[j])
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k, arg)
	c(0, 0)

#true value around (21.146, 16.027, 2.724, 0.574, 0.000)
z0 <- rexp(sum(dimx)+sum(dimlam))


funSSR(z0, dimx, dimlam, grobj=grobj, myarg, constr=g, myarg, grconstr=grg, myarg, compl=phiFB, echo=TRUE)

jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, myarg, grconstr=grg, myarg, 
heconstr=heg, myarg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, myarg, grconstr=grg, myarg, heconstr=heg, myarg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1, maxit=itermax))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, myarg, grconstr=grg, myarg, heconstr=heg, myarg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1, maxit=itermax))



#-------------------------------------------------------------------------------
# (4) Example of GNE with 4 solutions(!)
#-------------------------------------------------------------------------------

myarg <- list(C=c(2, 3), D=c(4,0))

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	dij <- 1*(i == j)
	other <- ifelse(i == 1, 2, 1)
	2*(x[i] - arg$C[i])*(x[other] - arg$D[i])^4*dij + 4*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^3*(1-dij) 
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
	dij <- 1*(i == j)
	dik <- 1*(i == k)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[other] - arg$D[i])^4*dij*dik + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*dij*(1-dik)
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*(1-dij)*dik
	res + 12*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^2*(1-dij)*(1-dik)
}

dimlam <- c(1, 1)
#g_i(x)
g <- function(x, i)
	ifelse(i == 1, sum(x[1:2]) - 1, 2*x[1]+x[2]-2)
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	ifelse(i == 1, 1, 1 + 1*(i == j))
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0



z0 <- rexp(sum(dimx)+sum(dimlam))

funSSR(z0, dimx, dimlam, grobj=grobj, myarg, constr=g, grconstr=grg, compl=phiFB, echo=FALSE)

jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB, echo=FALSE)


GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=1, maxit=itermax))

GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Broyden", 
	control=list(trace=1, maxit=itermax))


#-------------------------------------------------------------------------------
# (5) another Example
#-------------------------------------------------------------------------------

# associated objective functions

dimx <- c(2, 2, 3)
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



# constraint linked functions

dimlam <- c(1, 2, 2)
#constraint function g_i(x)
g <- function(x, i)
{
	x <- x[1:7]
	if(i == 1)
		res <- sum( x^(1:7) ) -7
	if(i == 2)
		res <- c(sum(x) + prod(x) - 14, 20 - sum(x))
	if(i == 3)
		res <- c(sum(x^2) + 1, 100 - sum(x))
	res
}


#Gr_x_j g_i(x)
grfullg <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
		grad <- (1:7) * x ^ (0:6)
	if(i == 2)
	{
		grad <- 1 + sapply(1:7, function(i) prod(x[-i]))
		grad <- cbind(grad, -1)
	}
	if(i == 3)
		grad <- cbind(2*x, -1)
	
	if(i == 1)
		res <- grad[j]	
	if(i != 1)
		res <- grad[j,]	
	as.numeric(res)
}

#Gr_x_k Gr_x_j g_i(x)
hefullg <- function(x, i, j, k)
{
	x <- x[1:7]
	if(i == 1)
		he1 <- diag( c(0, 2, 6, 12, 20, 30, 42) * x ^ c(0, 0, 1:5) )
	if(i == 2)
	{
		he1 <- matrix(0, 7, 7)
		he1[1, -1] <- sapply(2:7, function(i) prod(x[-c(1, i)]))
		he1[2, -2] <- sapply(c(1, 3:7), function(i) prod(x[-c(2, i)]))
		he1[3, -3] <- sapply(c(1:2, 4:7), function(i) prod(x[-c(3, i)]))
		he1[4, -4] <- sapply(c(1:3, 5:7), function(i) prod(x[-c(4, i)]))
		he1[5, -5] <- sapply(c(1:4, 6:7), function(i) prod(x[-c(5, i)]))
		he1[6, -6] <- sapply(c(1:5, 7:7), function(i) prod(x[-c(6, i)]))
		he1[7, -7] <- sapply(1:6, function(i) prod(x[-c(7, i)]))
		
		he2 <- matrix(0, 7, 7)
	}
	if(i == 3)
	{
		he1 <- diag(rep(2, 7))
		he2 <- matrix(0, 7, 7)
	}
	if(i != 1)
		return( c(he1[j, k], he2[j, k])	)
	else				
		return( he1[j, k] )
}




# (3) compute Phi
#

z <- rexp(sum(dimx) + sum(dimlam))

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]

funSSR(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, compl=phiFB)

jacSSR(z, dimx, dimlam, heobj=hefullob, constr=g, grconstr=grfullg, 
	heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)

x <- GNE.nseq(z, dimx, dimlam, grobj=grfullob, NULL, heobj=hefullob, NULL, 
	constr=g, NULL, grconstr=grfullg, NULL, heconstr=hefullg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=0, maxit=itermax))

print(x)
summary(x)



