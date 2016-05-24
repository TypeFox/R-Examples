# Chapter 1.3 R commands

# Section 1.3.2

x=rnorm(10,mean=50,sd=10)
y=rnorm(10,mean=50,sd=10)
m=length(x)
n=length(y)
sp=sqrt(((m-1)*sd(x)^2+(n-1)*sd(y)^2)/(m+n-2))
t.stat=(mean(x)-mean(y))/(sp*sqrt(1/m+1/n))

tstatistic=function(x,y)
{
m=length(x)
n=length(y)
sp=sqrt(((m-1)*sd(x)^2+(n-1)*sd(y)^2)/(m+n-2))
t.stat=(mean(x)-mean(y))/(sp*sqrt(1/m+1/n))
return(t.stat)
}

data.x=c(1,4,3,6,5)
data.y=c(5,4,7,6,10)
tstatistic(data.x, data.y)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.3.3

# simulation algorithm for normal populations

alpha=.1; m=10; n=10 # sets alpha, m, n
N=10000 # sets the number of simulations
n.reject=0 # counter of num. of rejections
for (i in 1:N)
{
x=rnorm(m,mean=0,sd=1) # simulates xs from population 1
y=rnorm(n,mean=0,sd=1) # simulates ys from population 2
t.stat=tstatistic(x,y) # computes the t statistic
if (abs(t.stat)>qt(1-alpha/2,n+m-2))
   n.reject=n.reject+1 # reject if |t| exceeds critical pt
}
true.sig.level=n.reject/N # est. is proportion of rejections

s=readline(prompt="Type  <Return>   to continue : ")

# simulation algorithm for normal and exponential populations
# storing the values of the t statistic in vector tstat

m=10; n=10 
my.tsimulation=function()
  tstatistic(rnorm(m,mean=10,sd=2), rexp(n,rate=1/10))
tstat.vector=replicate(10000, my.tsimulation())

plot(density(tstat.vector),xlim=c(-5,8),ylim=c(0,.4),lwd=3)
curve(dt(x,df=18),add=TRUE)
