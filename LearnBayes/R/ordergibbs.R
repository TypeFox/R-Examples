ordergibbs=function(data,m)
{
# implements Gibbs sampling for table of means
# with prior belief in order restriction
# input:  data = data matrix with two columns [sample mean, sample size]
#         m = number of iterations of Gibbs sampling
# output:  matrix of simulated values of means where each row
#          represents one simulated draw

#####################################################
rnormt=function(n,mu,sigma,lo,hi)
{
# simulates n random variates from a normal(mu,sigma)
# distribution truncated on the interval (lo, hi)

p=pnorm(c(lo,hi),mu,sigma)
return(mu+sigma*qnorm(runif(n)*(p[2]-p[1])+p[1]))
}
#####################################################

y=data[,1]  # sample means
n=data[,2]  # sample sizes
s=.65       # assumed value of sigma for this example
I=8; J=5    # number of rows and columns in matrix

# placing vectors y, n into matrices

y=t(array(y,c(J,I)))
n=t(array(n,c(J,I)))
y=y[seq(8,1,by=-1),]
n=n[seq(8,1,by=-1),]

# setting up the matrix of values of the population means mu
# two rows and two columns are added that help in the simulation
# of individual values of mu from truncated normal distributions

mu0=Inf*array(1,c(I+2,J+2))
mu0[1,]=-mu0[1,]
mu0[,1]=-mu0[,1]
mu0[1,1]=-mu0[1,1]
mu=mu0

# starting value of mu that satisfies order restriction

m1=c(2.64,3.02,3.02,3.07,3.34)
m2=c(2.37,2.63,2.74,2.76,2.91)
m3=c(2.37,2.47,2.64,2.66,2.66)
m4=c(2.31,2.33,2.33,2.33,2.33)
m5=c(2.04,2.11,2.11,2.33,2.33)
m6=c(1.85,1.85,1.85,2.10,2.10)
m7=c(1.85,1.85,1.85,1.88,1.88)
m8=c(1.59,1.59,1.59,1.67,1.88)
muint=rbind(m8,m7,m6,m5,m4,m3,m2,m1)
mu[2:(I+1),2:(J+1)]=muint

MU=array(0,c(m,I*J))  # arry MU stores simulated values of mu

#####################  main loop #######################
for (k in 1:m)
{
	for (i in 2:(I+1))
	{
		for (j in 2:(J+1))
		{
			lo=max(c(mu[i-1,j],mu[i,j-1]))
			hi=min(c(mu[i+1,j],mu[i,j+1]))
			mu[i,j]=rnormt(1,y[i-1,j-1],s/sqrt(n[i-1,j-1]),lo,hi)	
		}
	}
	mm=mu[2:(I+1),2:(J+1)]
	MU[k,]=array(mm,c(1,I*J))
}

return(MU)
}
