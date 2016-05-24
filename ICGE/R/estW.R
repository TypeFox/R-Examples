estW<-function(d, dx0, pert="onegroup"){
################# ESTIMATION OF INCA STATISTIC ################################
# Input:
# d: distance matrix between individuals (nxn) 
# dx0: vector of length n with distancies from the specific individual to the
#      individuals of different groups.
# pert: integer vector indicating the group each individual belongs to.
#
# Output: value of INCA statistic (W), and values of projections (U1, ..., Uk)
###############################################################################

d <- as.matrix(d)
n<-dim(d)[1]
if (pert[1]=="onegroup"){pert <- rep(1,n)}
pert <- as.integer(pert)
k<-max(pert)

# populations must be named with numbers from 1 to k
if (length(tabulate(as.factor(pert))) != k)
  stop("Partitions must be named by factors or with numbers from 1 to k.")
# 0 can not be a partitions name
if (any(pert==0))
  stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")


# We need the geometrical variabilities
var <- vgeo(d,pert) 

delta <- deltas(d, pert)

phi <- proxi(d, dx0, pert)



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
#library(MASS) (Dagoenekoz ezarrita dago DESCRIPTION-en Depends atalean)
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

if (W<0){
    W<-0
    cat("Warning: possibly a negative value of INCA statistic.")
}

for (i in 1:k){
U[i]=phi[i]-W
}
out=list(Wvalue=W,  Uvalue=U)

class(out) <- "incaest"

return(out)
}
