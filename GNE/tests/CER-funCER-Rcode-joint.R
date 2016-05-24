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



# (3) compute H
#

z <- rexp(sum(dimx) + 2*sum(dimmu))

n <- sum(dimx)
m <- 0
x <- z[1:n]
lam <- NULL
mu <- z[(n+m+1):(n+m+dimmu)]
w <- z[(n+m+dimmu+1):(n+2*m+2*dimmu)]

resphi <- funCER(z, dimx, dimlam, grobj=grfullob, joint=h, grjoint=grh, dimmu=dimmu, echo=TRUE)


check <- c(grfullob(x, 1, 1) + mu %*% grh(x, 1),
	grfullob(x, 1, 2) + mu %*% grh(x, 2),
	grfullob(x, 2, 3) + mu %*% grh(x, 3),
	grfullob(x, 2, 4) + mu %*% grh(x, 4),
	grfullob(x, 3, 5) + mu %*% grh(x, 5),
	grfullob(x, 3, 6) + mu %*% grh(x, 6),
	grfullob(x, 3, 7) + mu %*% grh(x, 7),
	h(x) + w, 
	mu * w)

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
	
resjacphi <- jacCER(z, dimx, dimlam, heobj=hefullob, 
	joint=h, grjoint=grh, hejoint=heh, dimmu=dimmu)

	
#check
cat("\n\n________________________________________\n\n")


cat("\n\npart A\n\n")	


cat("\n\npart A\n\n")	

lam <- rep(0, 5)

checkA <- 
rbind(
c(hefullob(x, 1, 1, 1) + mu%*%heh(x, 1, 1), 
hefullob(x, 1, 1, 2) + mu%*%heh(x, 1, 2),
hefullob(x, 1, 1, 3) + mu%*%heh(x, 1, 3),
hefullob(x, 1, 1, 4) + mu%*%heh(x, 1, 4),
hefullob(x, 1, 1, 5) + mu%*%heh(x, 1, 5),
hefullob(x, 1, 1, 6) + mu%*%heh(x, 1, 6),
hefullob(x, 1, 1, 7) + mu%*%heh(x, 1, 7)
),
c(hefullob(x, 1, 2, 1) + mu%*%heh(x, 2, 1), 
hefullob(x, 1, 2, 2) + mu%*%heh(x, 2, 2),
hefullob(x, 1, 2, 3) + mu%*%heh(x, 2, 3),
hefullob(x, 1, 2, 4) + mu%*%heh(x, 2, 4),
hefullob(x, 1, 2, 5) + mu%*%heh(x, 2, 5),
hefullob(x, 1, 2, 6) + mu%*%heh(x, 2, 6),
hefullob(x, 1, 2, 7) + mu%*%heh(x, 2, 7)
),
c(hefullob(x, 2, 3, 1) + mu%*%heh(x, 3, 1), 
hefullob(x, 2, 3, 2) + mu%*%heh(x, 3, 2),
hefullob(x, 2, 3, 3) + mu%*%heh(x, 3, 3),
hefullob(x, 2, 3, 4) + mu%*%heh(x, 3, 4),
hefullob(x, 2, 3, 5) + mu%*%heh(x, 3, 5),
hefullob(x, 2, 3, 6) + mu%*%heh(x, 3, 6),
hefullob(x, 2, 3, 7) + mu%*%heh(x, 3, 7)
),
c(hefullob(x, 2, 4, 1) + mu%*%heh(x, 4, 1), 
hefullob(x, 2, 4, 2) + mu%*%heh(x, 4, 2), 
hefullob(x, 2, 4, 3) + mu%*%heh(x, 4, 3), 
hefullob(x, 2, 4, 4) + mu%*%heh(x, 4, 4), 
hefullob(x, 2, 4, 5) + mu%*%heh(x, 4, 5), 
hefullob(x, 2, 4, 6) + mu%*%heh(x, 4, 6), 
hefullob(x, 2, 4, 7) + mu%*%heh(x, 4, 7)
),
c(hefullob(x, 3, 5, 1) + mu%*%heh(x, 5, 1),  
hefullob(x, 3, 5, 2) + mu%*%heh(x, 5, 2),  
hefullob(x, 3, 5, 3) + mu%*%heh(x, 5, 3),  
hefullob(x, 3, 5, 4) + mu%*%heh(x, 5, 4),  
hefullob(x, 3, 5, 5) + mu%*%heh(x, 5, 5),  
hefullob(x, 3, 5, 6) + mu%*%heh(x, 5, 6),  
hefullob(x, 3, 5, 7) + mu%*%heh(x, 5, 7)
),
c(hefullob(x, 3, 6, 1) + mu%*%heh(x, 6, 1),   
hefullob(x, 3, 6, 2) + mu%*%heh(x, 6, 2),  
hefullob(x, 3, 6, 3) + mu%*%heh(x, 6, 3),  
hefullob(x, 3, 6, 4) + mu%*%heh(x, 6, 4),  
hefullob(x, 3, 6, 5) + mu%*%heh(x, 6, 5),  
hefullob(x, 3, 6, 6) + mu%*%heh(x, 6, 6),  
hefullob(x, 3, 6, 7) + mu%*%heh(x, 6, 7)
),
c(hefullob(x, 3, 7, 1) + mu%*%heh(x, 7, 1),   
hefullob(x, 3, 7, 2) + mu%*%heh(x, 7, 2),  
hefullob(x, 3, 7, 3) + mu%*%heh(x, 7, 3),  
hefullob(x, 3, 7, 4) + mu%*%heh(x, 7, 4),  
hefullob(x, 3, 7, 5) + mu%*%heh(x, 7, 5),  
hefullob(x, 3, 7, 6) + mu%*%heh(x, 7, 6),  
hefullob(x, 3, 7, 7) + mu%*%heh(x, 7, 7)
)
)


print(resjacphi[1:n, 1:n] - checkA)


cat("\n\n________________________________________\n\n")
cat("\n\npart B\n\n")	


checkB <- rbind(grh(x, 1), grh(x, 2), grh(x, 3), grh(x, 4), grh(x, 5), grh(x, 6), grh(x, 7))

print(resjacphi[1:n, (n+1):(n+m+dimmu)] - checkB)	

cat("\n\n________________________________________\n\n")
cat("\n\npart C\n\n")	

checkC <- matrix(0, n, m+dimmu)

print(resjacphi[1:n, (n+m+dimmu+1):(n+2*m+2*dimmu)] - checkC)

cat("\n\n________________________________________\n\n")
cat("\n\npart D\n\n")	

checkD <- cbind(grh(x, 1), grh(x, 2), grh(x, 3), grh(x, 4), grh(x, 5), grh(x, 6), grh(x, 7))


print(resjacphi[(n+1):(n+m+dimmu), 1:n] - checkD)


cat("\n\n________________________________________\n\n")
cat("\n\npart E\n\n")	


checkE <- matrix(0, m+dimmu, m+dimmu)

print(resjacphi[(n+1):(n+m+dimmu), (n+1):(n+m+dimmu)] - checkE)


cat("\n\n________________________________________\n\n")
cat("\n\npart F\n\n")	


checkF <- diag(m+dimmu)

print(resjacphi[(n+1):(n+m+dimmu), (n+m+dimmu+1):(n+2*m+2*dimmu)] - checkF)


cat("\n\n________________________________________\n\n")
cat("\n\npart G\n\n")	


checkG <- matrix(0, m+dimmu, n)

print(resjacphi[(n+m+dimmu+1):(n+2*m+2*dimmu), 1:n] - checkG)

cat("\n\n________________________________________\n\n")
cat("\n\npart H\n\n")	


checkH <- diag(w)

print(resjacphi[(n+m+dimmu+1):(n+2*m+2*dimmu), (n+1):(n+m+dimmu)] - checkH)


cat("\n\n________________________________________\n\n")
cat("\n\npart I\n\n")	


checkI <- diag(mu)

print(resjacphi[(n+m+dimmu+1):(n+2*m+2*dimmu), (n+m+dimmu+1):(n+2*m+2*dimmu)] - checkI)

checkjack <- rbind(cbind(checkA, checkB, checkC), cbind(checkD, checkE, checkF), cbind(checkG, checkH, checkI))

if(sum(abs(checkjack - resjacphi)) > .Machine$double.eps^(2/3))
	stop("wrong result")
