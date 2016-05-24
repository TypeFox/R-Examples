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

dimlam <- dimw <- c(1, 2, 2)

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




# (3) compute H
#

z <- rexp(sum(dimx) + sum(dimlam) + sum(dimw))

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]
w <- z[(n+m+1):(n+m+m)]

resphi <- funCER(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, echo=TRUE)


check <- c(grfullob(x, 1, 1) + lam[1] * grfullg(x, 1, 1),
	grfullob(x, 1, 2) + lam[1] * grfullg(x, 1, 2),
	grfullob(x, 2, 3) + lam[2:3] %*% grfullg(x, 2, 3),
	grfullob(x, 2, 4) + lam[2:3] %*% grfullg(x, 2, 4),
	grfullob(x, 3, 5) + lam[4:5] %*% grfullg(x, 3, 5),
	grfullob(x, 3, 6) + lam[4:5] %*% grfullg(x, 3, 6),
	grfullob(x, 3, 7) + lam[4:5] %*% grfullg(x, 3, 7),
	g(x, 1) + w[1], 
	g(x, 2) + w[2:3], 
	g(x, 3) + w[4:5], 
	lam * w)

#check
cat("\n\n________________________________________\n\n")

#part A
print(cbind(check, res=as.numeric(resphi))[1:length(x), ])
#part B
print(cbind(check, res=as.numeric(resphi))[(length(x)+1):length(z), ])
	

if(sum(abs(check - resphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")



# (4) compute Jac H
#
	
resjacphi <- jacCER(z, dimx, dimlam, heobj=hefullob, constr=g, grconstr=grfullg, 
	heconstr=hefullg)

	
#check
cat("\n\n________________________________________\n\n")


cat("\n\npart A\n\n")	


checkA <- 
rbind(
c(hefullob(x, 1, 1, 1) + lam[1]*hefullg(x, 1, 1, 1), 
hefullob(x, 1, 1, 2) + lam[1]*hefullg(x, 1, 1, 2),
hefullob(x, 1, 1, 3) + lam[1]*hefullg(x, 1, 1, 3),
hefullob(x, 1, 1, 4) + lam[1]*hefullg(x, 1, 1, 4),
hefullob(x, 1, 1, 5) + lam[1]*hefullg(x, 1, 1, 5),
hefullob(x, 1, 1, 6) + lam[1]*hefullg(x, 1, 1, 6),
hefullob(x, 1, 1, 7) + lam[1]*hefullg(x, 1, 1, 7)
),
c(hefullob(x, 1, 2, 1) + lam[1]*hefullg(x, 1, 2, 1), 
hefullob(x, 1, 2, 2) + lam[1]*hefullg(x, 1, 2, 2),
hefullob(x, 1, 2, 3) + lam[1]*hefullg(x, 1, 2, 3),
hefullob(x, 1, 2, 4) + lam[1]*hefullg(x, 1, 2, 4),
hefullob(x, 1, 2, 5) + lam[1]*hefullg(x, 1, 2, 5),
hefullob(x, 1, 2, 6) + lam[1]*hefullg(x, 1, 2, 6),
hefullob(x, 1, 2, 7) + lam[1]*hefullg(x, 1, 2, 7)
),
c(hefullob(x, 2, 3, 1) + lam[2:3]%*%hefullg(x, 2, 3, 1), 
hefullob(x, 2, 3, 2) + lam[2:3]%*%hefullg(x, 2, 3, 2),
hefullob(x, 2, 3, 3) + lam[2:3]%*%hefullg(x, 2, 3, 3),
hefullob(x, 2, 3, 4) + lam[2:3]%*%hefullg(x, 2, 3, 4),
hefullob(x, 2, 3, 5) + lam[2:3]%*%hefullg(x, 2, 3, 5),
hefullob(x, 2, 3, 6) + lam[2:3]%*%hefullg(x, 2, 3, 6),
hefullob(x, 2, 3, 7) + lam[2:3]%*%hefullg(x, 2, 3, 7)
),
c(hefullob(x, 2, 4, 1) + lam[2:3]%*%hefullg(x, 2, 4, 1), 
hefullob(x, 2, 4, 2) + lam[2:3]%*%hefullg(x, 2, 4, 2), 
hefullob(x, 2, 4, 3) + lam[2:3]%*%hefullg(x, 2, 4, 3), 
hefullob(x, 2, 4, 4) + lam[2:3]%*%hefullg(x, 2, 4, 4), 
hefullob(x, 2, 4, 5) + lam[2:3]%*%hefullg(x, 2, 4, 5), 
hefullob(x, 2, 4, 6) + lam[2:3]%*%hefullg(x, 2, 4, 6), 
hefullob(x, 2, 4, 7) + lam[2:3]%*%hefullg(x, 2, 4, 7)
),
c(hefullob(x, 3, 5, 1) + lam[4:5]%*%hefullg(x, 3, 5, 1),  
hefullob(x, 3, 5, 2) + lam[4:5]%*%hefullg(x, 3, 5, 2),  
hefullob(x, 3, 5, 3) + lam[4:5]%*%hefullg(x, 3, 5, 3),  
hefullob(x, 3, 5, 4) + lam[4:5]%*%hefullg(x, 3, 5, 4),  
hefullob(x, 3, 5, 5) + lam[4:5]%*%hefullg(x, 3, 5, 5),  
hefullob(x, 3, 5, 6) + lam[4:5]%*%hefullg(x, 3, 5, 6),  
hefullob(x, 3, 5, 7) + lam[4:5]%*%hefullg(x, 3, 5, 7)
),
c(hefullob(x, 3, 6, 1) + lam[4:5]%*%hefullg(x, 3, 6, 1),   
hefullob(x, 3, 6, 2) + lam[4:5]%*%hefullg(x, 3, 6, 2),  
hefullob(x, 3, 6, 3) + lam[4:5]%*%hefullg(x, 3, 6, 3),  
hefullob(x, 3, 6, 4) + lam[4:5]%*%hefullg(x, 3, 6, 4),  
hefullob(x, 3, 6, 5) + lam[4:5]%*%hefullg(x, 3, 6, 5),  
hefullob(x, 3, 6, 6) + lam[4:5]%*%hefullg(x, 3, 6, 6),  
hefullob(x, 3, 6, 7) + lam[4:5]%*%hefullg(x, 3, 6, 7)
),
c(hefullob(x, 3, 7, 1) + lam[4:5]%*%hefullg(x, 3, 7, 1),   
hefullob(x, 3, 7, 2) + lam[4:5]%*%hefullg(x, 3, 7, 2),  
hefullob(x, 3, 7, 3) + lam[4:5]%*%hefullg(x, 3, 7, 3),  
hefullob(x, 3, 7, 4) + lam[4:5]%*%hefullg(x, 3, 7, 4),  
hefullob(x, 3, 7, 5) + lam[4:5]%*%hefullg(x, 3, 7, 5),  
hefullob(x, 3, 7, 6) + lam[4:5]%*%hefullg(x, 3, 7, 6),  
hefullob(x, 3, 7, 7) + lam[4:5]%*%hefullg(x, 3, 7, 7)
)
)


print(resjacphi[1:n, 1:n] - checkA)


cat("\n\n________________________________________\n\n")
cat("\n\npart B\n\n")	


checkB <- 
rbind(
cbind(c(grfullg(x, 1, 1), grfullg(x, 1, 2)), c(0, 0), c(0, 0), c(0, 0), c(0, 0)),
cbind(c(0, 0), rbind(grfullg(x, 2, 3), grfullg(x, 2, 4)), c(0, 0), c(0, 0)),
cbind(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), rbind(grfullg(x, 3, 5), grfullg(x, 3, 6), grfullg(x, 3, 7)))
)

print(resjacphi[1:n, (n+1):(n+m)] - checkB)	

cat("\n\n________________________________________\n\n")
cat("\n\npart C\n\n")	

checkC <- matrix(0, n, m)

print(resjacphi[1:n, (n+m+1):(n+2*m)] - checkC)

cat("\n\n________________________________________\n\n")
cat("\n\npart D\n\n")	




checkD <- 
 t(
cbind(
rbind(
grfullg(x, 1, 1) ,
grfullg(x, 1, 2) ,
grfullg(x, 1, 3) ,
grfullg(x, 1, 4) ,
grfullg(x, 1, 5) ,
grfullg(x, 1, 6) ,
grfullg(x, 1, 7) 
),
rbind(
grfullg(x, 2, 1) ,
grfullg(x, 2, 2) ,
grfullg(x, 2, 3) ,
grfullg(x, 2, 4) ,
grfullg(x, 2, 5) ,
grfullg(x, 2, 6) ,
grfullg(x, 2, 7) 
),
rbind(
grfullg(x, 3, 1) ,
grfullg(x, 3, 2) ,
grfullg(x, 3, 3) ,
grfullg(x, 3, 4) ,
grfullg(x, 3, 5) ,
grfullg(x, 3, 6) ,
grfullg(x, 3, 7) 
)
)
)



print(resjacphi[(n+1):(n+m), 1:n] - checkD)


cat("\n\n________________________________________\n\n")
cat("\n\npart E\n\n")	


checkE <- matrix(0, m, m)

print(resjacphi[(n+1):(n+m), (n+1):(n+m)] - checkE)


cat("\n\n________________________________________\n\n")
cat("\n\npart F\n\n")	


checkF <- diag(m)

print(resjacphi[(n+1):(n+m), (n+m+1):(n+2*m)] - checkF)


cat("\n\n________________________________________\n\n")
cat("\n\npart G\n\n")	


checkG <- matrix(0, m, n)

print(resjacphi[(n+m+1):(n+2*m), 1:n] - checkG)

cat("\n\n________________________________________\n\n")
cat("\n\npart H\n\n")	


checkH <- diag(w)

print(resjacphi[(n+m+1):(n+2*m), (n+1):(n+m)] - checkH)


cat("\n\n________________________________________\n\n")
cat("\n\npart I\n\n")	


checkI <- diag(lam)

print(resjacphi[(n+m+1):(n+2*m), (n+m+1):(n+2*m)] - checkI)

checkjack <- rbind(cbind(checkA, checkB, checkC), cbind(checkD, checkE, checkF), cbind(checkG, checkH, checkI))

if(sum(abs(checkjack - resjacphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")
