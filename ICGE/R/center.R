center <- function(dat){
############## Centers and standardizes the configurations  ##################
#      such that trace(X'i*Xi)=1,  before the procrustes analysis 
# Input:
# dat: array with configurations
# Output:
# dat: centered and standardized data.
##############################################################################

nrep<-dim(dat)[3]
tpoints<-dim(dat)[1]
p<-dim(dat)[2]
datini<-dat

if (is.na(nrep)){  #When no a third dimension, theres only one configuration
   nrep <- 1
   dat<- array(dat, dim=c(tpoints,p,1))
}


# Centering matrix
one=rep(1,tpoints)
J=one %*% t(one)
I=diag(tpoints) #identity matrix 
H=I-(1/tpoints)*J #centering matrix

# Centered configurations
for (i in 1:nrep){
dat[,,i]<- H %*% dat[,,i]
}


# Standarize data between different configurations 
# Trace (Xi' Xi)=trace(Xj' Xj)=1

for (i in 1:nrep)
{
A<-0
  for (j in 1:p){
     A <- A + var(dat[,j,i])
  }

a<-sqrt((tpoints-1)*A)
dat[,,i] <- dat[,,i]/a
}

if (is.na(dim(datini)[3])){  
dat<-dat[,,1]
}


return(dat)

}
