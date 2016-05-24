estW_simple<-function(var, delta, phi){
######################## Estimation of INCA W (faster version) ###############
# Input:
# var: vector of geometric variabilities
# delta: matrix of distances between clusters
# phi: proximity vector
# Output:
# Wvalue: estimation of W
# Uvalue: projections U1, ..., Uk
##############################################################################

k <- length(var)
U<-matrix(0,k,1)
M<- matrix(0,k-1,k-1)
N<- matrix(0,k-1,1)

for (i in 1:(k-1)){
    for (j in i:(k-1)){

        M[i,j]=delta[i,k]+delta[j,k]-delta[i,j]

    }
}

for (i in 1:(k-1)){
    for (j in i:(k-1)){
        M[j,i]=M[i,j]
    }
}

for (i in 1:(k-1)){
    N[i]=delta[i,k]+phi[k]-phi[i]
}
#library(MASS)
alpha<-ginv(M)%*%N

aux <- sum(alpha)
alpha <- c(alpha, 1-aux)


waux <- 0
for (i in 1:(k-1)){
    for(j in (i+1):k){
        waux <- waux + alpha[i]*alpha[j]*delta[i,j]
    }
}

W <- sum(alpha*phi) - waux


if ( W < 0 ){
    W <- 0
}

for (i in 1:k){
U[i]=phi[i]-W
}
out=list(Wvalue=W, Uvalue=U)
return(out)
}
