parameters <-
function(N,ns,a=0.1,E=0.99,p=0.99,method="default"){
if (missing(ns)) {
if (method == "default"){
ns=floor(N*0.5)+1
}
else if (method == "small.k"){
ns=floor(N*a)+1
}
else if (method == "appro.k"){
para2=function(ns,N,a=0.1,E=0.99,p=0.99){
ns=ceiling(ns)
n=N-floor(N*a)
r=ceiling(log(1-E)/(log(n-ns)-log(n)))
pg=choose(n,ns)/choose(N,ns)
f=function(l,r,pg){
t=0
for (i in 1:r){
t=t+choose(l,i-1)*pg^(i-1)*(1-pg)^(l+1-i)
}
return=1-t-p
return
}
fc=function(a,b){
return=choose(a,b)-6.601574e+274
return
}
if (fc(choose(N,ns),r) < 6.601574e+274){
maxk=choose(N,ns)
}
else {
maxk=uniroot(fc,c(1,choose(N,ns)),b=r)$root
}
k=ceiling(uniroot(f,c(1,maxk),r=r,pg=pg)$root)
return=r+k
return
}
ns=floor(optimize(para2,c(floor(N*a)+1,floor(N*0.5)+1),N=N,a=a,E=E,p=p)$minimum)
}
}
n=N-floor(N*a)
r=ceiling(log(1-E)/(log(n-ns)-log(n)))
pg=choose(n,ns)/choose(N,ns)
f=function(l,r,pg){
t=0
for (i in 1:r){
t=t+choose(l,i-1)*pg^(i-1)*(1-pg)^(l+1-i)
}
return=1-t-p
return
}
fc=function(a,b){
return=choose(a,b)-6.601574e+274
return
}
if (fc(choose(N,ns),r) < 6.601574e+274){
maxk=choose(N,ns)
}
else {
maxk=uniroot(fc,c(1,choose(N,ns)),b=r)$root
}
k=ceiling(uniroot(f,c(1,maxk),r=r,pg=pg)$root)
list (ns=ns,
      r=r,
      k=k)
}
