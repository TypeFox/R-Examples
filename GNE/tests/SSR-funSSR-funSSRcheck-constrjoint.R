if(!require("GNE"))stop("this test requires package GNE.")


# (1) associated objective functions
#

dimx <- c(2, 2, 3)

#Gr_x_j O_i(x)
grfullob <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- 3*(x - 1:7)^2
	}
	if(i == 2)
	{
		grad <- 1:7*(x - 1:7)^(0:6)
	}
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
	{
		he <- diag( 6*(x - 1:7) )
	}
	if(i == 2)
	{
		he <- diag( c(0, 2, 6, 12, 20, 30, 42)*(x - 1:7)^c(0, 0:5) )
	}
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



# (2) constraint linked functions
#

dimlam <- c(1, 2, 2)

#constraint function g_i(x)
g <- function(x, i)
{
	x <- x[1:7]
	#cat(x[1:5], "|", i, "\n")
	if(i == 1)
		res <- sum( x^(1:7) ) -7
	if(i == 2)
		res <- c(sum(x) + prod(x) - 14, 20 - sum(x))
	if(i == 3)
		res <- c(sum(x^2) + 1, 100 - sum(x))
	#cat("res", res + par$a, "\n")	
	res
}


#Gr_x_j g_i(x)
grfullg <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- (1:7) * x ^ (0:6)
	}
	if(i == 2)
	{
		grad <- 1 + sapply(1:7, function(i) prod(x[-i]))
		grad <- cbind(grad, -1)
	}
	if(i == 3)
	{
		grad <- cbind(2*x, -1)
	}
	if(i == 1) res <- grad[j]	
	if(i != 1) res <- grad[j,]	
	as.numeric(res)
}



#Gr_x_k Gr_x_j g_i(x)
hefullg <- function(x, i, j, k)
{
	x <- x[1:7]
	if(i == 1)
	{
		he1 <- diag( c(0, 2, 6, 12, 20, 30, 42) * x ^ c(0, 0, 1:5) )
	}
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




dimmu <- 3

#constraint function h(x)
h <- function(x)
{
	x <- x[1:7]
	c(prod(x) - 1, sum(x^2) -2, sum(x^3) -3)
}
grh <- function(x, j)
{
	x <- x[1:7]
	c(prod(x[-j]), 2*x[j], 3*x[j]^2)
}
heh <- function(x, j, k)
{
	x <- x[1:7] 
	c(prod(x[-c(j,k)]), 2*(j==k), 6*x[j]*(j==k))
}



# (3) compute Phi
#

z <- rexp(sum(dimx) + sum(dimlam) + dimmu)

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]
mu <- z[-(1:(n+m))]

resphi <- GNE:::funSSRcheck(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, compl=phiFB, joint=h, grjoint=grh, dimmu=dimmu)


check <- c(grfullob(x, 1, 1) + lam[1] * grfullg(x, 1, 1) + mu %*% grh(x, 1),
	grfullob(x, 1, 2) + lam[1] * grfullg(x, 1, 2) + mu %*% grh(x, 2),
	grfullob(x, 2, 3) + lam[2:3] %*% grfullg(x, 2, 3) + mu %*% grh(x, 3),
	grfullob(x, 2, 4) + lam[2:3] %*% grfullg(x, 2, 4) + mu %*% grh(x, 4),
	grfullob(x, 3, 5) + lam[4:5] %*% grfullg(x, 3, 5) + mu %*% grh(x, 5),
	grfullob(x, 3, 6) + lam[4:5] %*% grfullg(x, 3, 6) + mu %*% grh(x, 6),
	grfullob(x, 3, 7) + lam[4:5] %*% grfullg(x, 3, 7) + mu %*% grh(x, 7),
	phiFB(-g(x, 1), lam[1]), 
	phiFB( -g(x, 2)[1], lam[2]), 
	phiFB( -g(x, 2)[2], lam[3]), 
	phiFB( -g(x, 3)[1], lam[4]), 
	phiFB( -g(x, 3)[2], lam[5]),
	phiFB( -h(x)[1], mu[1]), 
	phiFB( -h(x)[2], mu[2]),
	phiFB( -h(x)[3], mu[3]))
	

#check
cat("\n\n________________________________________\n\n")

#part A
print(cbind(check, res=as.numeric(resphi))[1:n, ])
#part B
print(cbind(check, res=as.numeric(resphi))[(n+1):length(z), ])


if(sum(abs(check - resphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")



z <- rexp(sum(dimx) + dimmu)
n <- sum(dimx)
m <- 0
x <- z[1:n]
mu <- z[-(1:(n+m))]

resphi <- GNE:::funSSRcheck(z, dimx, grobj=grfullob, compl=phiFB, joint=h, grjoint=grh, dimmu=dimmu)


check <- c(grfullob(x, 1, 1) + mu %*% grh(x, 1),
	grfullob(x, 1, 2) + mu %*% grh(x, 2),
	grfullob(x, 2, 3) + mu %*% grh(x, 3),
	grfullob(x, 2, 4) + mu %*% grh(x, 4),
	grfullob(x, 3, 5) + mu %*% grh(x, 5),
	grfullob(x, 3, 6) + mu %*% grh(x, 6),
	grfullob(x, 3, 7) + mu %*% grh(x, 7),
	phiFB( -h(x)[1], mu[1]), 
	phiFB( -h(x)[2], mu[2]),
	phiFB( -h(x)[3], mu[3]))
	

#check
cat("\n\n________________________________________\n\n")

#part A
print(cbind(check, res=as.numeric(resphi))[1:n, ])
#part B
print(cbind(check, res=as.numeric(resphi))[(n+1):length(z), ])


if(sum(abs(check - resphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")







# (4) compute Jac Phi
#

z <- rexp(sum(dimx) + sum(dimlam) + dimmu)

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]
mu <- z[-(1:(n+m))]


resjacphi <- jacSSR(z, dimx, dimlam, heobj=hefullob, constr=g, 	
	grconstr=grfullg, heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB,
	joint=h, grjoint=grh, hejoint=heh, dimmu=dimmu)


	
resjaccheck <- GNE:::jacSSRcheck(z, dimx, dimlam, heobj=hefullob, constr=g, 	
	grconstr=grfullg, heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB,
	joint=h, grjoint=grh, hejoint=heh, dimmu=dimmu)

#check
cat("\n\n________________________________________\n\n")
cat("\n\npart A\n\n")	

print(resjacphi[1:n, 1:n] - resjaccheck[1:n, 1:n])


cat("\n\n________________________________________\n\n")
cat("\n\npart B\n\n")	

print(resjacphi[1:n, (n+1):(n+m)] - resjaccheck[1:n, (n+1):(n+m)])	



cat("\n\n________________________________________\n\n")
cat("\n\npart C\n\n")	

print(resjacphi[1:n, (n+m+1):(n+m+dimmu)] - resjaccheck[1:n, (n+m+1):(n+m+dimmu)])
	


cat("\n\n________________________________________\n\n")
cat("\n\npart D\n\n")	


print(resjacphi[(n+1):(n+m), 1:n] - resjaccheck[(n+1):(n+m), 1:n])


cat("\n\n________________________________________\n\n")
cat("\n\npart E\n\n")	
 

print(resjacphi[(n+1):(n+m), (n+1):(n+m)] - resjaccheck[(n+1):(n+m), (n+1):(n+m)])

cat("\n\n________________________________________\n\n")
cat("\n\npart F\n\n")	


print(resjacphi[(n+1):(n+m), (n+m+1):(n+m+dimmu)] - resjaccheck[(n+1):(n+m), (n+m+1):(n+m+dimmu)])

cat("\n\n________________________________________\n\n")
cat("\n\npart G\n\n")	

print(resjacphi[(n+m+1):(n+m+dimmu), 1:n] - resjaccheck[(n+m+1):(n+m+dimmu), 1:n])

cat("\n\n________________________________________\n\n")
cat("\n\npart H\n\n")	

print(resjacphi[(n+m+1):(n+m+dimmu), (n+1):(n+m)] - resjaccheck[(n+m+1):(n+m+dimmu), (n+1):(n+m)])


cat("\n\n________________________________________\n\n")
cat("\n\npart I\n\n")	

print(resjacphi[(n+m+1):(n+m+dimmu), (n+m+1):(n+m+dimmu)] - resjaccheck[(n+m+1):(n+m+dimmu), (n+m+1):(n+m+dimmu)])



if(sum(abs(resjacphi - resjaccheck)) > .Machine$double.eps^(2/3))
	stop("wrong result")









z <- rexp(sum(dimx) + dimmu)

resjacphi <- jacSSR(z, dimx, heobj=hefullob, 
	gcompla=GrAphiFB, gcomplb=GrBphiFB,
	joint=h, grjoint=grh, hejoint=heh, dimmu=dimmu)

resjaccheck <- GNE:::jacSSRcheck(z, dimx, heobj=hefullob, 
	gcompla=GrAphiFB, gcomplb=GrBphiFB,
	joint=h, grjoint=grh, hejoint=heh, dimmu=dimmu)

if(sum(abs(resjacphi - resjaccheck)) > .Machine$double.eps^(2/3))
	stop("wrong result")
	
	
z <- rexp(sum(dimx) + sum(dimlam))
	
resjacphi <- jacSSR(z, dimx, dimlam, heobj=hefullob, constr=g, 	
	grconstr=grfullg, heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)

resjaccheck <- GNE:::jacSSRcheck(z, dimx, dimlam, heobj=hefullob, constr=g, 	
	grconstr=grfullg, heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


if(sum(abs(resjacphi - resjaccheck)) > .Machine$double.eps^(2/3))
	stop("wrong result")

