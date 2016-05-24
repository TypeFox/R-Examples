reconstruct=function(niter,y)
{
numb=dim(y)[1]
x=0*y

mu=matrix(0,niter,6)
sigma2=rep(0,niter)

mu[1,]=c(35,50,65,84,92,120)
sigma2[1]=100

beta=rep(1,niter)
xcum=matrix(0,numb^2,6)
n=rep(0,6)


dali=c(6667.729,7245.159,7856.514,8523.00,9242.127,10025.211,10896.380,
11877.379,12985.344,14360.080,16062.470,18408.592,22755.124,33163.207,
35947.756,36745.675,38286.608,38534.912,38531.211,38916.662,38495.781)
thefunc=approxfun(seq(0,2,length=21),dali)

# the vector dali comes from a crude if time
# consumming approximation of the normalizing constant
# based on 21 points

for (i in 2:niter)
{

lvr=0

for (k in 1:numb)
{
for (l in 1:numb)
{
n[1]=xneig4(x,k,l,1)
n[2]=xneig4(x,k,l,2)
n[3]=xneig4(x,k,l,3)
n[4]=xneig4(x,k,l,4)
n[5]=xneig4(x,k,l,5)
n[6]=4-n[1]-n[2]-n[3]-n[4]-n[5]

x[k,l]=sample(1:6,1,prob=exp(beta[i-1]*n)*dnorm(y[k,l],mu[i-1,],sqrt(sigma2[i-1])))
xcum[(k-1)*100+l,x[k,l]]=xcum[(k-1)*100+l,x[k,l]]+1

lvr=lvr+n[x[k,l]]
}
}

mu[i,1]=truncnorm(1,mean(y[x==1]),sqrt(sigma2[i-1]/sum(x==1)),0,mu[i-1,2])
mu[i,2]=truncnorm(1,mean(y[x==2]),sqrt(sigma2[i-1]/sum(x==2)),mu[i,1],mu[i-1,3])
mu[i,3]=truncnorm(1,mean(y[x==3]),sqrt(sigma2[i-1]/sum(x==3)),mu[i,2],mu[i-1,4])
mu[i,4]=truncnorm(1,mean(y[x==4]),sqrt(sigma2[i-1]/sum(x==4)),mu[i,3],mu[i-1,5])
mu[i,5]=truncnorm(1,mean(y[x==5]),sqrt(sigma2[i-1]/sum(x==4)),mu[i,4],mu[i-1,6])
mu[i,6]=truncnorm(1,mean(y[x==6]),sqrt(sigma2[i-1]/sum(x==5)),mu[i,5],255)

sese=sum((y-mu[i,1])^2*(x==1)+(y-mu[i,2])^2*(x==2)+(y-mu[i,3])^2*(x==3)+(y-mu[i,4])^2*(x==4)+(y-mu[i,5])^2*(x==5)+(y-mu[i,6])^2*(x==6))
sigma2[i]=1/rgamma(1,numb^2/2,sese/2)

betatilde=beta[i-1]+runif(1,-0.05,0.05)

laccept=lvr*(betatilde-beta[i-1])+integrate(thefunc,betatilde,beta[i-1])$value
if (runif(1)<=exp(laccept)) beta[i]=betatilde else beta[i]=beta[i-1]

# print(round(c(i,beta[i],mu[i,],sigma2[i]),2))
}

list(beta=beta,mu=mu,sigma2=sigma2,xcum=xcum)
}
