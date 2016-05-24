whitenoise.test <-
function(x)
 {
 # This programme performs the Lobato-Velasco white noise test 
 # Econometric Theory, Vol. 20, Issue 04, 2004
 n<- length(x)
 # compute acf function
 gam<- acf(x,type="covariance")
 # sample autocovariance at lag 0 is gam0
 gam0<- gam$acf[1]
 # compute periodogram of x
 ILAM<- spec.pgram(x,taper=0,fast=FALSE)
 ILAM<- ILAM$spec
 T<- length(ILAM)
 # compute CVM statistic MN
 P2<- (T^(-1))*(sum(ILAM^2))
 MN<- (P2/gam0^2)-1
 # compute test statistic tMN
 tMN<- sqrt(T)*(MN-1)
 print("no. of observations")
 print(n)
 print("T")
 print(T)
 print("CVM stat MN")
 print(MN)
 print("tMN")
 print(tMN)
 pval<- pnorm(tMN,mean=0,sd=2,lower.tail=FALSE)
 # 2 tail test against N(0,4)
 # H0 the data are white noise
 test<- pval*2
 if (test > 1) test<- 2-test
 print("test value")
 print(test)
 }

