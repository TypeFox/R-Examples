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




# (3) compute Phi
#

z <- rexp(sum(dimx) + sum(dimlam))

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]

resphi <- funSSR(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, compl=phiFB)

check <- GNE:::funSSRcheck(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, compl=phiFB)



#check
cat("\n\n________________________________________\n\n")

#part A
print(cbind(check, res=as.numeric(resphi))[1:length(x), ])
#part B
print(cbind(check, res=as.numeric(resphi))[(length(x)+1):length(z), ])


if(sum(abs(check - resphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")


# (4) compute Jac Phi
#
	
resjacphi <- jacSSR(z, dimx, dimlam, heobj=hefullob, constr=g, grconstr=grfullg, 
	heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


checkjac <- GNE:::jacSSRcheck(z, dimx, dimlam, heobj=hefullob, constr=g, grconstr=grfullg, 
heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)


#check
cat("\n\n________________________________________\n\n")
cat("\n\npart A\n\n")	

print(resjacphi[1:n, 1:n] - checkjac[1:n, 1:n])


cat("\n\n________________________________________\n\n")
cat("\n\npart B\n\n")	

print(resjacphi[1:n, (n+1):(n+m)] - checkjac[1:n, (n+1):(n+m)])	


cat("\n\n________________________________________\n\n")
cat("\n\npart C\n\n")	


print(resjacphi[(n+1):(n+m), 1:n] - checkjac[(n+1):(n+m), 1:n])


cat("\n\n________________________________________\n\n")

cat("\n\npart D\n\n")	


print(resjacphi[(n+1):(n+m), (n+1):(n+m)] - checkjac[(n+1):(n+m), (n+1):(n+m)])


if(sum(abs(checkjac - resjacphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")
