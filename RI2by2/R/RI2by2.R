#Install packages
#install.packages('gtools',repos='http://cran.case.edu/')
#install.packages('compiler')

#Robins (1988)
Robins.CI = function(data,level) {
m0 = data[2,1]/(data[2,1]+data[2,2])
m1 = data[1,1]/(data[1,1]+data[1,2])
n = sum(data)
tau.hat = m1-m0
if (m1 >= m0) {
 r = ((2*m0-m1)*(1-m1)-m0*(1-m0))/n}
if (m1 < m0) {
 r = ((2*m1-m0)*(1-m0)-m1*(1-m1))/n}
se = sqrt(m0*(1-m0)/(data[2,1]+data[2,2])+m1*(1-m1)/(data[1,1]+data[1,2])+r)
lower = max(tau.hat-qnorm(1-level/2)*se,-1)
upper = min(tau.hat+qnorm(1-level/2)*se,1)
output.all = list(tau.hat=tau.hat,lower=lower,upper=upper)
return(output.all)
}

#Attributable effects approach
AE.CI = function(data,level) {
n = sum(data)
m = sum(data[1,])
tau.hat = data[1,1]/m-data[2,1]/(n-m)
  
f = double()
A1 = -data[1,2]:data[1,1]
for (i in (1:length(A1))) {
j = data
j[1,1]=j[1,1]-A1[i]
j[1,2]=j[1,2]+A1[i]
f[i] = fisher.test(j)$p.value
}
lower1 = min(A1[f>=(level/2)])
upper1 = max(A1[f>=(level/2)])

g = double()
A0 = -data[2,1]:data[2,2]
for (i in (1:length(A0))) {
j = data
j[2,1]=j[2,1]+A0[i]
j[2,2]=j[2,2]-A0[i]
g[i] = fisher.test(j)$p.value
}
lower2 = min(A0[g>=(level/2)])
upper2 = max(A0[g>=(level/2)])

#With Bonferroni adjustment
lower = (lower1+lower2)/n
upper = (upper1+upper2)/n

output.all = list(tau.hat=tau.hat,lower=lower,upper=upper)
return(output.all)
}

#Permutation of O(n^4)
#Exact re-randomization matrix for small n choose m
library('gtools')
nchoosem = function(n,m) {
c = choose(n,m)
trt = combinations(n,m)
Z = matrix(NA,c,n)
for (i in 1:c) {
Z[i,trt[i,]] = 1
Z[i,-trt[i,]] = 0
}
return(Z)
}

#Sample from re-randomization matrix
comb = function(n,m,nperm) {
trt = matrix(NA,nperm,m)
for (i in 1:nperm) {
trt[i,] = sample(n,m)
}
Z = matrix(NA,nperm,n)
for (i in 1:nperm) {
Z[i,trt[i,]] = 1
Z[i,-trt[i,]] = 0
}
return(Z)
}
#Permutation based p-value for each delta0 to test
#Corrected rounding issue 4/15/14
library('compiler')
pval2 = function(y.1,y.0,delta0,Z) {
 m = length(y.1[is.na(y.1)==0])
 n = m+length(y.0[is.na(y.0)==0])
 tau.hat = mean(y.1[is.na(y.1)==0])-mean(y.0[is.na(y.0)==0])
 dat = matrix(NA,n,3) 
 dat[,3] = delta0
 dat[1:m,1] = y.1[1:m]
 dat[(m+1):n,2] = y.0[(m+1):n]
 dat[1:m,2] = y.1[1:m]-delta0[1:m]
 dat[(m+1):n,1] = y.0[(m+1):n]+delta0[(m+1):n]  
 tau0 = mean(dat[,3])
 t.c = Z%*%dat[,1]/(m)-(1-Z)%*%dat[,2]/(n-m)
 p = mean(round(abs(t.c-tau0),15)>=round(abs(tau.hat-tau0),15))
 output.all = list(p,tau0,tau.hat,t.c)
 names(output.all) = list("p","tau0","tau.hat","t.c")
 return(output.all)
}
pval = cmpfun(pval2)

#Function to find permutation based confidence interval
Perm.CI2 = function(data,level,nperm) {
m = sum(data[1,])
n = m+sum(data[2,])
a = data[1,1]
b = data[1,2]
c = data[2,1]
d = data[2,2]
y.1 = c(rep(1,data[1,1]),rep(0,data[1,2]),rep(NA,n-m))
y.0 = c(rep(NA,m),rep(1,data[2,1]),rep(0,data[2,2]))
Z.obs = c(rep(1,m),rep(0,n-m))
Y.obs = c(rep(1,a),rep(0,b),grep(1,c),rep(0,d))
tau.hat = mean(y.1[is.na(y.1)==0])-mean(y.0[is.na(y.0)==0])

d1 = matrix(0,a,a)
d1[col(d1) >= row(d1)] = 1
d1 = cbind(rep(0,a),d1)

d2 = matrix(0,b,b)
d2[col(d2) >= row(d2)] = -1
d2 = cbind(rep(0,b),d2)

d3 = matrix(0,c,c)
d3[col(d3) >= row(d3)] = -1
d3 = cbind(rep(0,c),d3)

d4 = matrix(0,d,d)
d4[col(d4) >= row(d4)] = 1
d4 = cbind(rep(0,d),d4)

C = choose(n,m)
if (C<=nperm) Z = nchoosem(n,m) else Z = comb(n,m,nperm)

p = double()
tau0 = double()
for (g in 1:(a+1)) {
 for (h in 1:(b+1)) {
  for (i in 1:(c+1)) {
   for (j in 1:(d+1)) {
    k = pval(y.1,y.0,c(d1[,g],d2[,h],d3[,i],d4[,j]),Z)
    p = c(p,k$p)
    tau0 = c(tau0,k$tau0)
    #can also record delta here
   }
  }
 }
}

lower = min(tau0[p>=level])
upper = max(tau0[p>=level])
#output.all = list(tau.hat=tau.hat,lower=lower,upper=upper,tau0=tau0,p=p)
output.all = list(tau.hat=tau.hat,lower=lower,upper=upper)
#can also output p, tau0 here for conjecture
return(output.all)
}
Perm.CI = cmpfun(Perm.CI2)
