ci.prop <-
function(x,n,a){

v1=2*x
v2=2*(n-x+1)
v3=2*(x+1)
v4=2*(n-x)

a2=1-a/2
invb1=1
if(v1>0) invb1=qf(a2,v2,v1)

invb2=1
if(v4>0) invb2=qf(a2,v3,v4)
plexact=v1/(v1+v2*invb1)
puexact=v3*invb2/(v4+v3*invb2)

z=-qnorm(1-a2)
p=x/n;
plapprox=p-z*sqrt(p*(1-p)/n)
puapprox=p+z*sqrt(p*(1-p)/n)

result<-matrix(rep(0,4),ncol=2)
rownames(result)<-c("Exact CI","Approximate CI")
colnames(result)<-c("Lower Bound","Upper Bound")
result[1,1]<-plexact;result[1,2]<-puexact;result[2,1]<-plapprox;result[2,2]<-puapprox
result
}

