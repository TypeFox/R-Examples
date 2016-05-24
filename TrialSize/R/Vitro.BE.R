Vitro.BE <-
function(alpha, beta,delta, sigmaBT,sigmaBR,sigmaWT,sigmaWR,thetaBE){

U<-function(m,n,alpha, beta,delta, sigmaBT,sigmaBR,sigmaWT,sigmaWR,thetaBE)
{
U=((abs(delta)+qnorm(alpha)*sqrt(sigmaBT^2/m+sigmaBR^2/m))^2-delta^2)^2
+(sigmaBT^2+sigmaWT^2)^2*((m-1)/qchisq(1-alpha, m-1)-1)^2
+(1+thetaBE)^2*(sigmaBR^2+sigmaWR^2)^2*((m-1)/qchisq(alpha, m-1)-1)^2
}

sigmaT=sqrt(sigmaBT^2+sigmaWT^2)
sigmaR=sqrt(sigmaBR^2+sigmaWR^2)

gamma=delta^2+sigmaT^2-sigmaR^2-thetaBE*sigmaR^2

for (i in 1:1000){
bound=gamma+sqrt(U(i,1,alpha, 0.05,delta, sigmaBT,sigmaBR,sigmaWT,sigmaWR,thetaBE))
+sqrt(U(i,1,alpha, beta,delta, sigmaBT,sigmaBR,sigmaWT,sigmaWR,thetaBE))
print(c(i,bound))
}
}
