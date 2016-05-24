S_H1 <-
function(N,nr,nt,y,X,K,family,beta_new,Rhat,fihat) {

aindex=nr*cumsum(nt)
index=c(0,aindex[-length(aindex)])

eta=X%*%beta_new
mu=family$linkinv(eta)

sum201<-matrix(0,(K+1),1)      #gradient:S
sum301<-matrix(0,(K+1),(K+1))  #naive variance:H
sum401<-matrix(0,(K+1),(K+1))  #robust variance:M

for (i in 1:N) {
ym<-matrix(0,nt[i]*nr,1)
bigD<-matrix(0,nt[i]*nr,(K+1))
bigA<-matrix(0,nt[i]*nr,nt[i]*nr)
for (r in 1:nr)    {   
for (j in 1:nt[i]) {
#cat("r",r,"j",j,"\n")
ym[j+nt[i]*(r-1)]<- y[j+index[i]+nt[i]*(r-1)]-mu[j+index[i]+nt[i]*(r-1)] 
bigA[j+nt[i]*(r-1),j+nt[i]*(r-1)]<-family$variance(mu)[j+index[i]+nt[i]*(r-1)]
for (k in 1:(K+1)) {
bigD[j+nt[i]*(r-1),k]<-family$mu.eta(eta)[j+index[i]+nt[i]*(r-1)]*X[j+index[i]+nt[i]*(r-1),k]
#cat("i",i,"r",r,"j",j,"k",k,"\n")
} # for k
} # for r
} # for j

##working covariance matrix
bigV<-sqrt(bigA)%*%Rhat[1:(nr*nt[i]),1:(nr*nt[i]),i]%*%sqrt(bigA)
bigV<-fihat*bigV

##This is S in Wang et al.(2012)
sum200<-t(bigD)%*%ginv(bigV)%*%ym      #this is like gradient
sum201<-sum201+sum200

##This is H in Wang et al.(2012) and B in Inan (2015)
sum300<-t(bigD)%*%ginv(bigV)%*%bigD    #this is for information matrix.
sum301<-sum301+sum300

##This is B in Inan (2015)##
sum400<-t(bigD)%*%ginv(bigV)%*%ym%*%t(ym)%*%ginv(bigV)%*%bigD
sum401<-sum401+sum400

#cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
}  #end of i

S<-sum201
H<-sum301
M<-sum401

return(list("S"=S,"H"=H,"M"=M))
}
