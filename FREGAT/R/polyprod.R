# Function from package 'fda' (c) 2014
#  --------------------------------------------------------------

polyprod <- function(Coeff1, Coeff2){
# POLYCONV computes products of polynomials defined by columns of 
#   coefficient matrices Coeff1 and Coeff2

#  Last modified 6 June 2005

polyorder1 <- dim(Coeff1)[1] 
norder1    <- dim(Coeff1)[2] 
polyorder2 <- dim(Coeff2)[1] 
norder2    <- dim(Coeff2)[2] 
ndegree1 <- polyorder1 - 1 
ndegree2 <- polyorder2 - 1 

#  if the degrees are not equal, pad out the smaller matrix with 0s

if(ndegree1 != ndegree2){
    if(ndegree1 > ndegree2)  Coeff2 <- rbind(Coeff2, matrix(0,ndegree1-ndegree2,norder2))
    else                     Coeff1 <- rbind(Coeff1, matrix(0,ndegree2-ndegree1,norder1))
}

#  find order of the product
D <- max(c(ndegree1,ndegree2))   # maximum degree
N <- 2*D+1                       # order of product

#  compute the coefficients for the products
convmat <- array(0,c(norder1,norder2,N)) 
for (i in 0:(D-1)){
    ind <- c(0:i) + 1 
    if (length(ind) == 1) {
        convmat[,,i+1] <-     outer(Coeff1[ind,    ],Coeff2[i-ind+2,]) 
        convmat[,,N-i] <-     outer(Coeff1[D-ind+2,],Coeff2[D-i+ind,])	
    } else {
        convmat[,,i+1] <- crossprod(Coeff1[ind,    ],Coeff2[i-ind+2,]) 
        convmat[,,N-i] <- crossprod(Coeff1[D-ind+2,],Coeff2[D-i+ind,])
    }
}
ind <- c(0:D) + 1 
convmat[,,D+1] <-     crossprod(Coeff1[ind,    ],Coeff2[D-ind+2,])

if (ndegree1 != ndegree2) {
	convmat <- convmat[,,1:(ndegree1+ndegree2+1)] 
	convmat <- array(convmat,c(norder1,norder2,ndegree1+ndegree2+1))
}

return(convmat)
}

