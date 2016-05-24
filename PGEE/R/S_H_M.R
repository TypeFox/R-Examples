S_H_M<-function(N,nt,y,X,K,family,beta_new,Rhat,fihat) {

aindex=cumsum(nt)
index=c(0,aindex[-length(aindex)])

eta=X%*%beta_new
mu=family$linkinv(eta)

sum201<-matrix(0,(K+1),1)      #gradient:S
sum301<-matrix(0,(K+1),(K+1))  #naive variance:H
sum401<-matrix(0,(K+1),(K+1))  #robust variance:M

for (i in 1:N) {
ym<-matrix(0,nt[i],1)
bigD<-matrix(0,nt[i],(K+1))
bigA<-matrix(0,nt[i],nt[i])
for (j in 1:nt[i]) {
#cat("j",j,"\n")
ym[j]<- y[j+index[i]]-mu[j+index[i]] 
bigA[j,j]<-family$variance(mu)[j+index[i]]
for (k in 1:(K+1)) {
bigD[j,k]<-family$mu.eta(eta)[j+index[i]]*X[j+index[i],k]
#cat("i",i,"j",j,"k",k,"\n")
} # for k
} # for j

##working covariance matrix
bigV<-sqrt(bigA)%*%Rhat[1:nt[i],1:nt[i],i]%*%sqrt(bigA)
#like bigSigma
##bigV<-fihat*bigV

##This is S in Wang et al.(2012)
sum200<-t(bigD)%*%ginv(bigV)%*%ym      #this is like gradient
sum201<-sum201+sum200

##This is H in Wang et al.(2012) and A in Inan (2015)
sum300<-t(bigD)%*%ginv(bigV)%*%bigD    #this is for information matrix.
sum301<-sum301+sum300

##This is B in Inan (2015)##
sum400<-t(bigD)%*%ginv(bigV)%*%ym%*%t(ym)%*%ginv(bigV)%*%bigD
sum401<-sum401+sum400

#cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
}  #end of i

S<-fihat*sum201
H<-fihat*sum301
M<-fihat*sum401

return(list("S"=S,"H"=H,"M"=M))
}
