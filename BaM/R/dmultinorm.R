# Description: 	dmultinorm function, see page 393
# Usage:	dmultinorm (xval,yval,mu.vector,sigma.matrix)



dmultinorm <- function(xval,yval,mu.vector,sigma.matrix)  {
   normalizer <- (2*pi*sigma.matrix[1,1]*sigma.matrix[2,2]
                  *sqrt(1-sigma.matrix[1,2]^2))^(-1)
   like <- exp(-(1/(2*(1-sigma.matrix[1,2]^2)))* (
               ((xval-mu.vector[1])/sigma.matrix[1,1])^2
               -2*sigma.matrix[1,2]*(((xval-mu.vector[1])
               /sigma.matrix[1,1])*
               ((yval-mu.vector[2])/sigma.matrix[2,2]))
               +((yval-mu.vector[2])/sigma.matrix[2,2])^2 ))
   normalizer*like
}

