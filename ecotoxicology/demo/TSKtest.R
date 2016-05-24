# TSKTest
# Test cases for the Trimmed Spearman-Karber method, as per Hamilton and EPA.
# Written by Brenton R. Stone, last revised May 18 2010
# Translated to R by Jose Gama 2013

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
Trims<-function(x,p,n)
{

conf=erf(2/sqrt(2))

tResult=TSK(x,p,n,0,conf)
mu0<-tResult$mu
gsd<-tResult$gsd
low0<-tResult$left
up0<-tResult$right
tResult=TSK(x,p,n,0.05,conf)
mu5<-tResult$mu
gsd<-tResult$gsd
low5<-tResult$left
up5<-tResult$right
tResult=TSK(x,p,n,0.1,conf)
mu1<-tResult$mu
gsd<-tResult$gsd
low1<-tResult$left
up1<-tResult$right
tResult=TSK(x,p,n,0.2,conf)
mu2<-tResult$mu
gsd<-tResult$gsd
low2<-tResult$left
up2<-tResult$right

cat(sprintf('mid:\t%4.2f \t%4.2f \t%4.2f \t%4.2f\n',mu0,mu5,mu1,mu2))
cat(sprintf('\tlow:\t%4.2f \t%4.2f \t%4.2f \t%4.2f\n',low0,low5,low1,low2))
cat(sprintf('\thigh:\t%4.2f \t%4.2f \t%4.2f \t%4.2f\n',up0,up5,up1,up2))
}

#Hamilton's test cases.
#Table I
x1=c(15.54,20.47,27.92,35.98,55.52)
n1=c(20,20,20,19,20)
p1a=c(0,0,0,5.26,100)/100*n1
p1b=c(0,5,0,5.26,100)/100*n1
p1c=c(0,5,0,15.79,100)/100*n1
p1d=c(0,5,5,94.74,100)/100*n1
p1e=c(0,5,5,100,100)/100*n1

#Table V
x4=c(7.8,13,22,36,60,100)
n4=10
p4a=c(0,0,10,100,100,100)/100*n4
p4b=c(0,0,70,100,100,100)/100*n4
p4c=c(0,0,10,40,100,100)/100*n4
p4d=c(0,0,20,70,100,100)/100*n4
p4e=c(0,0,20,30,100,100)/100*n4

# uncomment the sink commands to save all the output in a text file
#sink('TSKtestResult.txt')
cat('Spearman-Karber, pct\n')
cat('\t\t0 \t5 \t10 \t20\n')
cat('1A\t')
Trims(x1,p1a,n1)
cat('\n1B\t')
Trims(x1,p1b,n1)
cat('\n1C\t')
Trims(x1,p1c,n1)
cat('\n1D\t')
Trims(x1,p1d,n1)
cat('\n1E\t')
Trims(x1,p1e,n1)
cat('\n4A\t')
Trims(x4,p4a,n4)
cat('\n4B\t')
Trims(x4,p4b,n4)
cat('\n4C\t')
Trims(x4,p4c,n4)
cat('\n4D\t')
Trims(x4,p4d,n4)
cat('\n4E\t')
Trims(x4,p4e,n4)

#From the documentation to the EPA's program for doing the TSK Method. See
#http://www.epa.gov/eerd/stat2.htm
cat('\nEPA:\t mu\t down\t up\n')
tResult=TSK(c(0,6.25,12.5,25,50,100),c(0,0,1,0,0,16),20,0.2,erf(2/sqrt(2)))
mu<-tResult$mu
gsd<-tResult$gsd
left<-tResult$left
right<-tResult$right
cat(sprintf('\t%4.2f\t%4.2f\t%4.2f\n',mu,left,right))

#sink()

