Propensity.Score.nostrata <-
function(alpha, beta,J,a,b,p1,phi){
p2=phi*p1/(phi*p1+1-p1)

p1star=sum(a*b*p1)/sum(a*b)
p2star=sum(a*(1-b)*p2)/sum(a*(1-b))
pstar=sum(a*(b*p1+(1-b)*p2))

b1=sum(a*b)
b2=sum(a*(1-b))

sigma1star=p1star*(1-p1star)/b1+p2star*(1-p2star)/b2
sigma0star=pstar*(1-pstar)*(1/b1+1/b2)

n=(sqrt(sigma0star)*qnorm(1-alpha/2)+sqrt(sigma1star)*qnorm(1-beta))^2/(p1star-p2star)^2
}
