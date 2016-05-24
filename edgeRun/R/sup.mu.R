sup.mu = function(phi,n){

newdata = data.frame(x=n*phi/2)
newdata[newdata>20] = 20

mu = predict(fit.2,newdata=newdata) #target phi_n=nphi_2/2

poisson.zone = newdata<0.2
mu[poisson.zone] = 22.98

m = (45e3 - 22.98) / (0.32 - 0.2)
C = -0.2*m + 22.98
linear.zone = newdata>=0.2 & newdata <0.32
mu[linear.zone] = m*newdata[linear.zone] + C

return(2*mu/n)
}