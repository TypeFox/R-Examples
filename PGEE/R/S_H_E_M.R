S_H_E_M <-
function(N,nt,y,X,K,family,beta_new,Rhat,fihat,lambda,pindex,eps=10^-6) {

aindex=cumsum(nt)
index=c(0,aindex[-length(aindex)])

eta=X%*%beta_new
mu=family$linkinv(eta)

#This is E on Wang et al.(2012)
E1<-diag(q_scad(abs(as.vector(beta_new)),lambda)/(abs(as.vector(beta_new))+eps))

if(is.null(pindex)==TRUE) {
E<-E1 
} else 
if(is.null(pindex)!=TRUE) {
E1[,pindex]<-0
E<-E1
}

sum201<-matrix(0,(K+1),1)      #gradient:S
sum301<-matrix(0,(K+1),(K+1))  #naive variance:H
sum401<-matrix(0,(K+1),(K+1))  #a component for robust variance:M

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
#bigV<-fihat*bigV

##This is S in Wang et al.(2012)
sum200<-t(bigD)%*%ginv(bigV)%*%ym      #this is like gradient
sum201<-sum201+sum200

##This is H in Wang et al.(2012)
sum300<-t(bigD)%*%ginv(bigV)%*%bigD    #this is for information matrix.
sum301<-sum301+sum300

##Speed up the code##
SSA=sqrt(ginv(bigA))
SRhat=ginv(Rhat[1:nt[i],1:nt[i],i])
SSAym=(SSA%*%ym)

sum400<-t(bigD)%*%SSA%*%SRhat%*%(SSAym%*%t(SSAym))%*%SRhat%*%SSA%*%bigD
sum401<-sum401+sum400

#cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
}  #end of i

S<-fihat*sum201
H<-fihat*sum301
E<-E
M<-fihat*sum401

return(list("S"=S,"H"=H,"E"=E,"M"=M))
}
