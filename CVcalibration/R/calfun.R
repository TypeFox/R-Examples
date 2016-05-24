calfun=function(x, y, CVx, CVy=CVx, lambda0=1)
{

n=length(x)

A1=cbind(c(1, mean(x)), c(mean(x), mean(x^2/(1+CVx^2))))
beta1=solve(A1)%*%c(mean(y), mean(x*y))
s1=rbind(as.vector(y-beta1[1]-beta1[2]*x), x*as.vector(y-beta1[1]-beta1[2]*x/(1+CVx^2)))
B1=s1%*%t(s1)/n
v1=solve(A1)%*%B1%*%solve(A1)/n
sigma1=sqrt(diag(v1))


A2=cbind(c(1, mean(y)), c(mean(x), mean(x*y)))
beta2=solve(A2)%*%c(mean(y), mean(y^2/(1+CVy^2)))
s2=rbind(as.vector(y-beta2[1]-beta2[2]*x), y*as.vector(y/(1+CVy^2)-beta2[1]-beta2[2]*x))
B2=s2%*%t(s2)/n
v2=solve(A2)%*%B2%*%solve(A2)/n
sigma2=sqrt(diag(v2))

s3=rbind(as.vector(y-beta2[1]-beta2[2]*x), x*as.vector(y-beta2[1]-beta2[2]*x/(1+CVx^2)), y*as.vector(y/(1+CVy^2)-beta2[1]-beta2[2]*x))
sigma3=solve(s3%*%t(s3))

objq=function(beta, sigma)
  {score=rbind(as.vector(y-beta[1]-beta[2]*x), x*as.vector(y-beta[1]-beta[2]*x/(1+CVx^2)), y*as.vector(y/(1+CVy^2)-beta[1]-beta[2]*x))
   score=apply(score,1,sum)
   result=t(score)%*%sigma%*%score
   return(result)
   }
beta3=optim(beta2, objq, sigma=sigma3)$par
s3=rbind(as.vector(y-beta3[1]-beta3[2]*x), x*as.vector(y-beta3[1]-beta3[2]*x/(1+CVx^2)), y*as.vector(y/(1+CVy^2)-beta3[1]-beta3[2]*x))

sigma3=solve(s3%*%t(s3))
beta3=optim(beta2, objq, sigma=sigma3)$par


s3=rbind(as.vector(y-beta3[1]-beta3[2]*x), x*as.vector(y-beta3[1]-beta3[2]*x/(1+CVx^2)), y*as.vector(y/(1+CVy^2)-beta3[1]-beta3[2]*x))

B3=s3%*%t(s3)/n
A3=cbind(c(1, mean(x)), c(mean(x), mean(x^2/(1+CVx^2))), c(mean(y), mean(x*y)))
v3=solve(A3%*%solve(B3)%*%t(A3))/n
sigma3=sqrt(diag(v3))


beta4=rep(0, 2)
a=cov(x, y)*(lambda0*var(x)+mean(x)^2)
b=(cov(x, y)^2-mean(x)^2*var(y))-lambda0*(cov(x, y)^2-var(x)*mean(y)^2)
c=-cov(x, y)*(var(y)+lambda0*mean(y)^2)
beta4[2]=(-b+sqrt(b^2-4*a*c))/2/a
beta4[1]=mean(y)-beta4[2]*mean(x)
rx2=beta4[2]*mean(x^2)/(mean(x*y)-beta4[1]*mean(x))
ry2=mean(y^2)/(beta4[1]*mean(y)+beta4[2]*mean(x*y))

A4=cbind(c(1, rx2*mean(x),  ry2*mean(y)), c(mean(x), mean(x^2), ry2*mean(x*y)), c(0, beta4[1]*mean(x)-mean(x*y), lambda0*beta4[1]*mean(y)+lambda0*beta4[2]*mean(x*y) ))
s4=rbind(as.vector(y-beta4[1]-beta4[2]*x), x*as.vector(y*rx2-beta4[1]*rx2-beta4[2]*x), y*as.vector(y-beta4[1]*ry2-beta4[2]*x*ry2))
B4=s4%*%t(s4)/n
v4=(solve(A4)%*%B4%*%t(solve(A4))/n)[1:2, 1:2]
sigma4=sqrt(diag(v4))

betax=beta1
betay=beta2
betaxy=beta3
betaratio=beta4
sex=sigma1
sey=sigma2
sexy=sigma3
seratio=sigma4

result=matrix(0, 4, 8)
result[1,]=c(betax, betay, betaxy, betaratio)
result[2,]=c(sex, sey, sexy, seratio)
result[3,]=c(betax-1.96*sex,  betay-1.96*sey, betaxy-1.96*sexy, betaratio-1.96*seratio)
result[4,]=c(betax+1.96*sex,  betay+1.96*sey, betaxy+1.96*sexy, betaratio+1.96*seratio)
colnames(result)=c("intercept-CVx only", "slope-CVx only", "intercept-CVy only", "slope-CVy only", "intercept-CVx and CVy", "slope-CVx and CVy", "intercept-CVy/CVx only", "slope-CVy/CVx")
rownames(result)=c("coef", "se", "lower CI", "upper CI")


return(result)

}
