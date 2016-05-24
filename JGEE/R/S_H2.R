S_H2 <-
function(N,nr,nt,y,X,K,family,beta_new,Rhat,fihat) {

aindex=nr*cumsum(nt)
index=c(0,aindex[-length(aindex)])

xnewmat=matrix(0,sum(nt)*nr,(nr*(K+1)))

for(i in 1:N){
xnewmat[(1+index[i]):(nt[i]*nr+index[i]),]=kronecker(diag(nr),X[(1+(i-1)*nt[i]):(i*nt[i]),])
#cat((1+index[i]):(nt[i]*nr+index[i]),"\n")
}

#X=xnewmat
eta=xnewmat%*%beta_new
#cat("beta_new",beta_new,"eta",eta,"\n")

if (is.character(family)) family <- get(family)
if (is.function(family))  family <- family()

mu=sapply(1:nr, function(r) family[[r]]$linkinv(eta[as.vector(sapply(1:N, function(i) sapply(1:nt[i], function(j) j+index[i]+nt[i]*(r-1))))]))

muy=matrix(0,sum(nt)*nr,1)
yy=matrix(0,sum(nt)*nr,1)

for(i in 1:N) {
for(j in 1:nt[i]){
for(r in 1:nr){
yy[j+index[i]+nt[i]*(r-1)]=y[j+(i-1)*nt[i],r]
muy[j+index[i]+nt[i]*(r-1)]=mu[j+(i-1)*nt[i],r]
}
}
}

mu=muy

#cat("yy", yy, "\n")
#cat("mu", mu, "\n")

sum201<-matrix(0,(nr*(K+1)),1)           #gradient
sum301<-matrix(0,(nr*(K+1)),(nr*(K+1)))  #naive variance
sum401<-matrix(0,(nr*(K+1)),(nr*(K+1)))  #robust variance

for (i in 1:N) {
ym<-matrix(0,nt[i]*nr,1)
bigD<-matrix(0,nt[i]*nr,(nr*(K+1)))
bigA<-matrix(0,nt[i]*nr,nt[i]*nr)
for (r in 1:nr)    {   
for (j in 1:nt[i]) {
#cat("r",r,"j",j,"\n")
ym[j+nt[i]*(r-1)]<- yy[j+index[i]+nt[i]*(r-1)]-mu[j+index[i]+nt[i]*(r-1)] 
bigA[j+nt[i]*(r-1),j+nt[i]*(r-1)]<-family[[r]]$variance(mu)[j+index[i]+nt[i]*(r-1)]
for (k in 1:(nr*(K+1))) {
bigD[j+nt[i]*(r-1),k]<- family[[r]]$mu.eta(eta)[j+index[i]+nt[i]*(r-1)]*xnewmat[j+index[i]+nt[i]*(r-1),k]
#cat("i",i,"r",r,"j",j,"k",k,"\n")
} # for k
} # for r
} # for j

#cat("i", i, "bigA", bigA, "\n")
#cat("i", i, "ym", ym, "\n")

##working covariance matrix
bigV<-sqrt(bigA)%*%Rhat[1:(nr*nt[i]),1:(nr*nt[i]),i]%*%sqrt(bigA)
bigV<-fihat*bigV

#cat("i", i, "bigV", bigV, "\n")

##This is S in Wang et al.(2012)
sum200<-t(bigD)%*%ginv(bigV)%*%ym      #this is like gradient
sum201<-sum201+sum200

##This is H in Wang et al.(2012) and B in Inan (2015)
sum300<-t(bigD)%*%ginv(bigV)%*%bigD    #this is for information matrix.
sum301<-sum301+sum300

##This is B in Inan (2015)##
sum400<-t(bigD)%*%ginv(bigV)%*%ym%*%t(ym)%*%ginv(bigV)%*%bigD
sum401<-sum401+sum400

#cat("i", i, "sum201", sum201, "\n")
#cat("i", i, "sum301", sum301, "\n")
#cat("i", i, "sum401", sum401, "\n")

#cat("i",i,"sum201",sum201,"sum301",sum301,"sum401",sum401,"\n")
}  #end of i

S<-sum201
H<-sum301
M<-sum401
Xx<-xnewmat

#cat("beta_new",beta_new,"S",S,"\n")

return(list("S"=S,"H"=H,"M"=M,"Xx"=xnewmat))
}
