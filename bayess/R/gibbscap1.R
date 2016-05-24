gibbscap1=function(nsimu,n1,c2,c3,N0=n1/runif(1),r10,r20){
# GIBBS SAMPLING FOR THE 2-STAGE OPEN POPULATION MODEL

N=rep(0,nsimu)
p=rep(0,nsimu)
q=rep(0,nsimu)
r1=rep(0,nsimu)
r2=rep(0,nsimu)
N[1]=N0
r1[1]=r10
r2[1]=r20
nplus=n1+c2+c3

for (i in 2:nsimu){

 uplus=N[i-1]-r1[i-1]-c2+n1-r1[i-1]-r2[i-1]-c3
 p[i]=rbeta(1,nplus+1,uplus+1)
 q[i]=rbeta(1,r1[i-1]+r2[i-1]+1,2*n1-2*r1[i-1]-r2[i-1]+1)
 N[i]=n1+rnbinom(1,n1,p[i])
 mm=min(n1-r2[i-1]-c3,n1-c2)
 pq=q[i]/(1+(1-q[i])*(1-p[i]))
 pr=lchoose(n1-c2,0:mm)+(0:mm)*log(pq)+lchoose(n1-(0:mm),r2[i-1]+c3)
 r1[i]=sample(0:mm,1,prob=exp(pr-max(pr)))
 r2[i]=rbinom(1,n1-r1[i]-c3,q[i]/(1+(1-q[i])*(1-p[i])))
 }

list(N=N,p=p,q=q,r1=r1,r2=r2)
}
