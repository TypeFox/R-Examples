BMDE <-
function(x,p)
{
k=length(x)
n_r=10000
for(m in 1:k)
{
if(x[m]<0)
{cat('Warning: arguments must be non-negative integers')
}
}

po=x+p
dr=rdirichlet(n_r,po)
a=0
l=0
u=0
dif=0 
for(j in 1:k)
{
a[j]=round(mean(dr[,j]),4)
l[j]=round(quantile(dr[,j],0.025),4)
u[j]=round(quantile(dr[,j],0.975),4)
dif[j]=u[j]-l[j]
}
cat('Mean\n')
print(a)
#
#
cat('Lower Limit\n')
print(l)
#
#
cat('Upper Limit\n')
print(u)
#
cat('Volume \n')
prod(dif)
}
