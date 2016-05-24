IBE <-
function(alpha, beta,delta, sigmaD,sigmaWT,sigmaWR,a,b,thetaIBE){
Sigma<-function(sigmaD,sigmaWT,sigmaWR,a,b)
{
Sigma=sigmaD^2+a*sigmaWT^2+b*sigmaWR^2
}


U<-function(n,alpha, beta,delta, sigmaD,sigmaWT,sigmaWR,a,b,thetaIBE)
{
U=((abs(delta)+qt(alpha,2*n-2)*Sigma(sigmaD,sigmaWT,sigmaWR,0.5,0.5)*sqrt(2/n)/2)^2-delta^2)^2
+Sigma(sigmaD,sigmaWT,sigmaWR,0.5,0.5)^4*((2*n-2)/qchisq(1-alpha, 2*n-2)-1)^2
+0.5^2*sigmaWT^4*((2*n-2)/qchisq(1-alpha, 2*n-2)-1)^2
+(1.5+thetaIBE)^2*sigmaWR^4*((2*n-2)/qchisq(alpha, 2*n-2)-1)^2
}

gamma=delta^2+sigmaD^2+sigmaWT^2-sigmaWR^2-thetaIBE*sigmaWR^2
for (i in 1:1000){
bound=gamma+sqrt(U(i,alpha, 0.05,delta, sigmaD,sigmaWT,sigmaWR,a,b,thetaIBE))
+sqrt(U(i,alpha, beta,delta, sigmaD,sigmaWT,sigmaWR,a,b,thetaIBE))
print(c(i,bound))
}
}
