
samplesize=function(x0, d0, x=seq(20, 200, length=1000), CVx, CVy)
{
lambda0=CVy^2/CVx^2

mu1=mean(x)
mu2=mean(x^2)/(1+CVx^2)
mu3=mean(x^3)/(1+3*CVx^2)
mu4=mean(x^4)/(1+6*CVx^2+3*CVx^4)

v11=(CVx^2+CVy^2)*mu2
v12=(CVy^2+(CVx^2-CVx^4)/(1+CVx^2))*mu3
v13=(CVx^2+(CVy^2-CVy^4)/(1+CVy^2))*mu3
v22=(CVy^2*(1+CVx^2)+CVx^2*(1+CVx^4)/(1+CVx^2)^2)*mu4
v23=((CVx^2-CVx^4)/(1+CVx^2)+(CVy^2-CVy^4)/(1+CVy^2)+CVx^2*CVy^2/(1+CVx^2)/(1+CVy^2))*mu4
v33=(CVx^2*(1+CVy^2)+CVy^2*(1+CVy^4)/(1+CVy^2)^2)*mu4


A1=cbind(c(1, mu1), c(mu1, mu2))
A2=cbind(c(1, mu1), c(mu1, mu2))
A3=cbind(c(1, mu1), c(mu1, mu2), c(mu1, mu2))
A4=cbind(c(1, (1+CVx^2)*mu1, (1+CVy^2)*mu1), c(mu1, (1+CVx^2)*mu2, (1+CVy^2)*mu2), c(0, -mu2, lambda0*mu2))
B1=cbind(c(v11, v12), c(v12, v22))
B2=cbind(c(v11, v13), c(v13, v33))
B3=cbind(c(v11, v12, v13), c(v12, v22, v23), c(v13, v23, v33))
B4=diag(c(1, 1+CVx^2, 1+CVy^2), 3, 3)%*%B3%*%diag(c(1, 1+CVx^2, 1+CVy^2), 3, 3)


sigma1=solve(A1)%*%B1%*%solve(A1)
sigma2=solve(A2)%*%B2%*%solve(A2)
sigma3=solve(A3%*%solve(B3)%*%t(A3))
sigma4=(solve(A4)%*%B4%*%t(solve(A4)))[1:2, 1:2]

x0=c(1, x0)
v1=t(x0)%*%sigma1%*%x0
v2=t(x0)%*%sigma2%*%x0
v3=t(x0)%*%sigma3%*%x0
v4=t(x0)%*%sigma4%*%x0


n1=(2*1.96)^2/d0^2*v1
n2=(2*1.96)^2/d0^2*v2
n3=(2*1.96)^2/d0^2*v3
n4=(2*1.96)^2/d0^2*v4

size=t(matrix(c(ceiling(n1), ceiling(n2), ceiling(n3), ceiling(n4))))
colnames(size)=c("CV of X only", "CV of Y only", "CV of both X and Y", "CVy/CVx")

return(list(size=size))
}

