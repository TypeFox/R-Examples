binegbin.glm=function(y1,y2,x1,x2){
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# author: Masakazu Iwasaki, Clinical Statistics Group, Schering-Plough K.K., Tokyo
# 
n=length(y1)
p1=ncol(x1)
p2=ncol(x2)
p=p1+p2
print(p)
eta1=eta(y1)
eta2=eta(y2)

beta1=qr.solve(t(x1)%*%x1)%*%t(x1)%*%eta1
beta2=qr.solve(t(x2)%*%x2)%*%t(x2)%*%eta2
beta=rbind(beta1,beta2)
print(beta)
mu1=inveta(x1%*%beta1)
mu2=inveta(x2%*%beta2)

# initial value
phi=sum((y1-mu1)^2+(y2-mu2)^2-(mu1+mu2))/sum(mu1^2+mu2^2)
alpha=sum((y1-mu1)*(y2-mu2))/(phi*sum(mu1*mu2))

for(j in 1:100){
FI=matrix(0,p,p)
Q=matrix(0,p,1)

for(i in 1:n){
xi1=c(x1[i,],rep(0,p2))
xi2=c(rep(0,p1),x2[i,])
xi=cbind(xi1,xi2)
yi=c(y1[i],y2[i])
mui=c(mu1[i],mu2[i])

#phii=((y1[i]-mu1[i])^2+(y2[i]-mu2[i])^2-(mu1[i]+mu2[i]))/(mu1[i]^2+mu2[i]^2)
#alphai=c(y1[i]-mu1[i])*(y2[i]-mu2[i])/(phii*mu1[i]*mu2[i])

fi=detadmu(mui)
FI=FI+xi%*%qr.solve(fi%*%sigma(mui,phi,alpha)%*%fi,tol=1e-20)%*%t(xi)
Q=Q+xi%*%qr.solve(sigma(mui,phi,alpha)%*%fi,tol=1e-20)%*%(yi-mui)
}

dbeta=qr.solve(FI,tol=1e-20)%*%Q
beta=beta+dbeta
beta1=beta[1:p1]
beta2=beta[(p1+1):(p1+p2)]
mu1=inveta(x1%*%beta1)
mu2=inveta(x2%*%beta2)

phi=sum((y1-mu1)^2+(y2-mu2)^2-(mu1+mu2))/sum(mu1^2+mu2^2)
alpha=sum((y1-mu1)*(y2-mu2))/(phi*sum(mu1*mu2))

lammda=(1+alpha)/(phi*alpha)
#print(phi)
#print(alpha)
#print(lammda)
}
#print(y1)
#print(mu1)
#print(y2)
#print(mu2)

#cormu=cor(mu1,mu2)
#cory=cor(y1,y2)
#print(cormu)
#print(cory)

dy1=y1-mu1
dy2=y2-mu2

#meany1_mean(y1)
#meany2_mean(y2)
meand1=mean(dy1)
meand2=mean(dy2)
stdd1=sd(dy1)
stdd2=sd(dy2)
meanmu1=mean(mu1)
meanmu2=mean(mu2)

mse1=sum(dy1^2)/(n-1)
mse2=sum(dy2^2)/(n-1)

#print(meany1)
#print(meany2)
#print(meand1)
#print(meand2)
#print(stdd1)
#print(stdd2)


asymcov=solve(FI)
list(beta1,beta2,phi,alpha,lammda,asymcov,meand1,meand2,mse1,mse2,meanmu1,meanmu2)
#list(beta1,beta2,phi,alpha,lammda,mu1,mu2,asymcov,meand1,meand2,mse1,mse2,meanmu1,meanmu2)
}

eta=function(mu){
log(mu)
}

inveta=function(eta){
exp(eta)
}

detadmu=function(mu){
diag(c(1/mu[1],1/mu[2]))
}

sigma=function(mu,phi,alpha){
matrix(c(mu[1]+phi*mu[1]^2,alpha*phi*mu[1]*mu[2],alpha*phi*mu[1]*mu[2],mu[2]+phi*mu[2]^2),2,2)
}