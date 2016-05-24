toCart <- function(a,b,c,alpha,beta,gamma, latticecoords) {
  alpha <- alpha * (pi/180)
  beta <- beta * (pi/180)
  gamma <- gamma * (pi/180)
  A <- matrix(ncol=3,nrow=3)
  v <- sqrt(1-cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma))
            
  A[2,1] <- A[3,1] <- A[3,2] <- 0
  A[1,1] <- a
  A[1,2] <- b*cos(gamma)
  A[2,2] <- b*sin(gamma)
  A[1,3] <- c*cos(beta)
  A[2,3] <- c * ((cos(alpha) - cos(beta) * cos(gamma) )/sin(gamma))
  A[3,3] <- c * (v/sin(gamma))
  res <- solve(A,t(latticecoords))
  res
            

}


base <- matrix(c(c(1/3,2/3,0),
                 c(2/3,1/3,1/2),
                 c(1/3,2/3,.37),
                 c(2/3,1/3,.87)
                 ), ncol=3, byrow=TRUE)
				 
toCart(4.2985,4.2985,7.0152,90,90,120, base)