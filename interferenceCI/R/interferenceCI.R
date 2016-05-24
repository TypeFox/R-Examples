#Load dependencies
library('gtools')
library('compiler')

#Example data set
#d = c(1,1,5,3,0,6,3,1,0,4,3,3,0,5,3,2,1,1,5,3,2,2,4,2,1,5,2,2,2,3,4,1,1,1,5,3,1,5,2,2)
#data.ex = array(d,c(2,2,10))
#assign.ex = c(1,0,0,0,1,1,0,1,1,0)

############
#Estimation#
############

######################################
#Estimates and bounds for the effects#
######################################
estbound = cmpfun(
function(g,data,m.a0,m.a1) {
k = length(g)
h = sum(g)
n = apply(data,3,sum)
m = apply(data,3,function(x) x[1,1]+x[1,2])
c11 = apply(data,3,function(x) x[1,1])
c12 = apply(data,3,function(x) x[1,2])
c21 = apply(data,3,function(x) x[2,1])
c22 = apply(data,3,function(x) x[2,2])
y0.a0 = double()
y1.a0 = double()
y0.a1 = double()
y1.a1 = double()
S = double()
Z = double()
Y = double()
w = double()
w.n = double()
w.c = double()
w.t = double()
w.a0 = double()
w.a1 = double()
for (i in 1:k) {
  S = c(S,rep(g[i],n[i]))
  Z = c(Z,c(rep(1,m[i]),rep(0,n[i]-m[i])))
  Y = c(Y,c(rep(1,c11[i]),rep(0,c12[i]),rep(1,c21[i]),rep(0,c22[i])))
  w = c(w,c(rep(1/m[i],m[i]),rep(1/(n[i]-m[i]),n[i]-m[i])))
  w.n = c(w.n,rep(1/n[i],n[i]))
  w.c = c(w.c,rep(1/(n[i]-m[i]),n[i]))
  w.t = c(w.t,rep(1/m[i],n[i]))
  w.a0 = c(w.a0,rep(m.a0[i]/n[i],n[i]))
  w.a1 = c(w.a1,rep(m.a1[i]/n[i],n[i]))
  y0.a0 = c(y0.a0,g[i]*rep(999,n[i])+(1-g[i])*(c(rep(999,m[i]),rep(1,c21[i]),rep(0,c22[i]))))
  y1.a0 = c(y1.a0,g[i]*rep(999,n[i])+(1-g[i])*(c(rep(1,c11[i]),rep(0,c12[i]),rep(999,n[i]-m[i]))))
  y0.a1 = c(y0.a1,(1-g[i])*rep(999,n[i])+g[i]*(c(rep(999,m[i]),rep(1,c21[i]),rep(0,c22[i]))))
  y1.a1 = c(y1.a1,(1-g[i])*rep(999,n[i])+g[i]*(c(rep(1,c11[i]),rep(0,c12[i]),rep(999,n[i]-m[i]))))  
}
y0.a0[y0.a0==999] = NA
y1.a0[y1.a0==999] = NA
y0.a1[y0.a1==999] = NA
y1.a1[y1.a1==999] = NA

#Compute estimators
y0.a0.hat = (1/(k-h))*sum(w*(1-S)*(1-Z)*Y)
y1.a0.hat = (1/(k-h))*sum(w*(1-S)*Z*Y)
y0.a1.hat = (1/h)*sum(w*S*(1-Z)*Y)
y1.a1.hat = (1/h)*sum(w*S*Z*Y)
y.a0.hat = (1/(k-h))*sum(w.n*(1-S)*Y)
y.a1.hat = (1/h)*sum(w.n*S*Y)
                       
DEa0.hat = y0.a0.hat-y1.a0.hat
DEa1.hat = y0.a1.hat-y1.a1.hat
IE.hat = y0.a0.hat-y0.a1.hat
TE.hat = y0.a0.hat-y1.a1.hat
OE.hat = y.a0.hat-y.a1.hat
                       
#Compute bounds
y0.a0.low = y0.a0
y0.a0.low[is.na(y0.a0.low)==1] = 0
y0.a0.high = y0.a0
y0.a0.high[is.na(y0.a0.high)==1] = 1
y0.a0.lb = (1/k)*sum(w.n*y0.a0.low)
y0.a0.ub = (1/k)*sum(w.n*y0.a0.high)

y1.a0.low = y1.a0
y1.a0.low[is.na(y1.a0.low)==1] = 0
y1.a0.high = y1.a0
y1.a0.high[is.na(y1.a0.high)==1] = 1
y1.a0.lb = (1/k)*sum(w.n*y1.a0.low)
y1.a0.ub = (1/k)*sum(w.n*y1.a0.high)

y0.a1.low = y0.a1
y0.a1.low[is.na(y0.a1.low)==1] = 0
y0.a1.high = y0.a1
y0.a1.high[is.na(y0.a1.high)==1] = 1
y0.a1.lb = (1/k)*sum(w.n*y0.a1.low)
y0.a1.ub = (1/k)*sum(w.n*y0.a1.high)                       

y1.a1.low = y1.a1
y1.a1.low[is.na(y1.a1.low)==1] = 0
y1.a1.high = y1.a1
y1.a1.high[is.na(y1.a1.high)==1] = 1
y1.a1.lb = (1/k)*sum(w.n*y1.a1.low)
y1.a1.ub = (1/k)*sum(w.n*y1.a1.high)

#Also for average effects
y.a0.lb = (1/k)*sum(w.n*(w.a0*y1.a0.low+(1-w.a0)*y0.a0.low))
y.a0.ub = (1/k)*sum(w.n*(w.a0*y1.a0.high+(1-w.a0)*y0.a0.high))

y.a1.lb = (1/k)*sum(w.n*(w.a1*y1.a1.low+(1-w.a1)*y0.a1.low))
y.a1.ub = (1/k)*sum(w.n*(w.a1*y1.a1.high+(1-w.a1)*y0.a1.high))

DEa0.lb = y0.a0.lb-y1.a0.ub
DEa0.ub = y0.a0.ub-y1.a0.lb

DEa1.lb = y0.a1.lb-y1.a1.ub
DEa1.ub = y0.a1.ub-y1.a1.lb

IE.lb = y0.a0.lb-y0.a1.ub
IE.ub = y0.a0.ub-y0.a1.lb

TE.lb = y0.a0.lb-y1.a1.ub
TE.ub = y0.a0.ub-y1.a1.lb

OE.lb = y.a0.lb-y.a1.ub
OE.ub = y.a0.ub-y.a1.lb

tab.eff = cbind(c(DEa0.lb,DEa1.lb,IE.lb,TE.lb,OE.lb),c(DEa0.hat,DEa1.hat,IE.hat,TE.hat,OE.hat),c(DEa0.ub,DEa1.ub,IE.ub,TE.ub,OE.ub))
colnames(tab.eff) = c("lower","est","upper")
rownames(tab.eff) = c("DEa0","DEa1","IE","TE","OE")

#tab.PO = cbind(c(y0.a0.lb,y1.a0.lb,y0.a1.lb,y1.a1.lb),c(y0.a0.hat,y1.a0.hat,y0.a1.hat,y1.a1.hat),c(y0.a0.ub,y1.a0.ub,y0.a1.ub,y1.a1.ub))
#colnames(tab.PO) = c("lower","est","upper")
#rownames(tab.PO) = c("y0.a0","y1.a0","y0.a1","y1.a1")

output.all = list(DEa0.hat=DEa0.hat,DEa1.hat=DEa1.hat,IE.hat=IE.hat,TE.hat=TE.hat,OE.hat=OE.hat,DEa0.lb=DEa0.lb,DEa0.ub=DEa0.ub,DEa1.lb=DEa1.lb,DEa1.ub=DEa1.ub,IE.lb=IE.lb,IE.ub=IE.ub,TE.lb=TE.lb,TE.ub=TE.ub,OE.lb=OE.lb,OE.ub=OE.ub,w.c=w.c,w.t=w.t,y0.a0=y0.a0,y1.a0=y1.a0,y0.a1=y0.a1,y1.a1=y1.a1,y0.a0.hat=y0.a0.hat,y1.a0.hat=y1.a0.hat,y0.a1.hat=y0.a1.hat,y1.a1.hat=y1.a1.hat,y0.a0.low=y0.a0.low,y0.a0.high=y0.a0.high,y1.a0.low=y1.a0.low,y1.a0.high=y1.a0.high,y0.a1.low=y0.a1.low,y0.a1.high=y0.a1.high,y1.a1.low=y1.a1.low,y1.a1.high=y1.a1.high,y0.a0.lb=y0.a0.lb,y0.a0.ub=y0.a0.ub,y1.a0.lb=y1.a0.lb,y1.a0.ub=y1.a0.ub,y0.a1.lb=y0.a1.lb,y0.a1.ub=y0.a1.ub,y1.a1.lb=y1.a1.lb,y1.a1.ub=y1.a1.ub,y.a0.lb=y.a0.lb,y.a0.ub=y.a0.ub,y.a1.lb=y.a1.lb,y.a1.ub=y.a1.ub,n=n,m=m,S=S,Z=Z,Y=Y,g=g,tab.eff=tab.eff)

return(output.all)
}
)

#estbound(assign.ex,data.ex,rep(3,10),rep(6,10))

###########
#Inference#
###########

##################################
#Hudgens Halloran (2008) approach#
##################################
HH = cmpfun(
function(eff,g,data,m.a0,m.a1,level) {
a = estbound(g,data,m.a0,m.a1)
n = a$n
h = sum(g)
k = length(g)

group = double()
for (i in 1:k) {
  group = c(group,rep(i,n[i]))
}

y0.a0.m = double()
y1.a0.m = double()
y.a0.m = double()
y0.a1.m = double()
y1.a1.m = double()
y.a1.m = double()  
y0.a0.v = double()
y1.a0.v = double()
y.a0.v = double()
y0.a1.v = double()
y1.a1.v = double()
y.a1.v = double()  
for (i in 1:k) {
 y0.a0.m = c(y0.a0.m,mean(a$Y[group==i & a$S==0 & a$Z==0]))
 y1.a0.m = c(y1.a0.m,mean(a$Y[group==i & a$S==0 & a$Z==1]))
 y.a0.m = c(y.a0.m,mean(a$Y[group==i & a$S==0]))
 y0.a1.m = c(y0.a1.m,mean(a$Y[group==i & a$S==1 & a$Z==0]))
 y1.a1.m = c(y1.a1.m,mean(a$Y[group==i & a$S==1 & a$Z==1]))
 y.a1.m = c(y.a1.m,mean(a$Y[group==i & a$S==1]))
 y0.a0.v = c(y0.a0.v,var(a$Y[group==i & a$S==0 & a$Z==0]))
 y1.a0.v = c(y1.a0.v,var(a$Y[group==i & a$S==0 & a$Z==1]))
 y.a0.v = c(y.a0.v,var(a$Y[group==i & a$S==0]))
 y0.a1.v = c(y0.a1.v,var(a$Y[group==i & a$S==1 & a$Z==0]))
 y1.a1.v = c(y1.a1.v,var(a$Y[group==i & a$S==1 & a$Z==1]))
 y.a1.v = c(y.a1.v,var(a$Y[group==i & a$S==1])) 
} 

DEa0.m = y0.a0.m-y1.a0.m
DEa1.m = y0.a1.m-y1.a1.m

effect = double()
error = double()

if (eff=='DEa0') {
 effect=a$DEa0.hat
 error = (h/(k*(k-h)))*var(DEa0.m,na.rm=TRUE)+(1/(k*(k-h)))*(sum((1/a$m)*y1.a0.v,na.rm=TRUE)+sum((1/(a$n-a$m))*y0.a0.v,na.rm=TRUE))
}

if (eff=='DEa1') {
 effect=a$DEa1.hat
 error = ((k-h)/(k*h))*var(DEa1.m,na.rm=TRUE)+(1/(k*h))*(sum((1/a$m)*y1.a1.v,na.rm=TRUE)+sum((1/(a$n-a$m))*y0.a1.v,na.rm=TRUE))
}

if (eff=='IE') {
 effect=a$IE.hat
 error = (1/(k-h))*var(y0.a0.m,na.rm=TRUE)+(1/h)*var(y0.a1.m,na.rm=TRUE)
}

if (eff=='TE') {
 effect=a$TE.hat
 error = (1/(k-h))*var(y0.a0.m,na.rm=TRUE)+(1/h)*var(y1.a1.m,na.rm=TRUE)
}

if (eff=='OE') {
 effect=a$OE.hat
 error = (1/(k-h))*var(y.a0.m,na.rm=TRUE)+(1/h)*var(y.a1.m,na.rm=TRUE)
}
 
lower.w = max(-1,effect-qnorm(1-level/2)*sqrt(error))
upper.w = min(1,effect+qnorm(1-level/2)*sqrt(error))
lower.ch = max(-1,effect-sqrt(error/level))
upper.ch = min(1,effect+sqrt(error/level))

output.all = list(est=effect,v=error,lower.w=lower.w,upper.w=upper.w,lower.ch=lower.ch,upper.ch=upper.ch)
return(output.all)
}
)
#Verify using Table 3 in HH08
#hh = array(c(16,18,12541-16,12541-18,26,54,11513-26,11513-54,17,119,10772-17,25134-119,22,122,8883-22,20727-122,15,92,5627-15,13130-92),c(2,2,5))
#e1 = HH('OE',c(1,1,0,0,0),hh,round(0.3*c(25082,23026,35906,29610,18757),0),round(0.5*c(25082,23026,35906,29610,18757),0),0.05)
#round(1000*e1$est,3)
#round(1000000*e1$v,3)

#Check for made up example k
#HH('OE',assign.ex,data.ex,rep(3,10),rep(6,10),0.05)


#############################
#Tchetgen VanderWeele (2012)#
#############################

TV = cmpfun(
function(eff,g,data,m.a0,m.a1,level) {
a = estbound(g,data,m.a0,m.a1)
n = a$n
q.0 = sum(1-g)/length(g)
q.1 = sum(g)/length(g)
N = length(g)
effect = double()
error = double()
if (eff=='DEa0') {
 effect=a$DEa0.hat
 l = double()
 s = double()
 for (i in 1:N) {
   l[i] = 2*(1-1/(choose(n[i],m.a0[i])))
   s[i] = (l[i]/q.0)^2
 }
 error = sqrt((4*(1/q.0-1)^2+sum(s)/N)*log(2/level)/(2*N))
}

if (eff=='DEa1') {
 effect=a$DEa1.hat
 l = double()
 s = double()
 for (i in 1:N) {
   l[i] = 2*(1-1/(choose(n[i],m.a1[i])))
   s[i] = (l[i]/q.1)^2
 }
 error = sqrt((4*(1/q.1-1)^2+sum(s)/N)*log(2/level)/(2*N))
}

if (eff=='IE') {
 effect=a$IE.hat
 l = double()
 s = double()
 for (i in 1:N) {
   l[i] = max((1/(q.0^2))*(1-1/choose(n[i],m.a0[i])),(1/(q.1^2))*(1-1/choose(n[i],m.a1[i])))
   s[i] = l[i]^2
 }
 error = sqrt((max(1/(q.0^2),1/(q.1^2))+sum(s)/N)*log(2/level)/(2*N))
}

if (eff=='TE') {
 effect=a$TE.hat
 l = double()
 s = double()
 for (i in 1:N) {
   l[i] = max((1/(q.0^2))*(1-1/choose(n[i],m.a0[i])),(1/(q.1^2))*(1-1/choose(n[i],m.a1[i])))
   s[i] = l[i]^2
 }
 error = sqrt((max(1/(q.0^2),1/(q.1^2))+sum(s)/N)*log(2/level)/(2*N))
}

if (eff=='OE') {
 effect=a$OE.hat
 l = double()
 s = double()
 for (i in 1:N) {
   l[i] = max((1/(q.0^2))*(1-1/choose(n[i],m.a0[i])),(1/(q.1^2))*(1-1/choose(n[i],m.a1[i])))
   s[i] = l[i]^2
 }
 error = sqrt((max(1/(q.0^2),1/(q.1^2))+sum(s)/N)*log(2/level)/(2*N))
}

lower = max(-1,effect-error)
upper = min(1,effect+error)
output.all = list(est=effect,v=error,lower=lower,upper=upper)
return(output.all)  
}
)

#Check for made up example k
#TV('OE',assign.ex,data.ex,rep(3,10),rep(6,10),0.05)

#Check for HH data
#TV('OE',c(1,1,0,0,0),hh,round(0.3*c(25082,23026,35906,29610,18757),0),round(0.5*c(25082,23026,35906,29610,18757),0),0.05)

#######################################################################
#Rigdon Hudgens (2014) exact inference using inverted permutation test#
#######################################################################

#generate targeted sharp null hypothesis
sample.n = cmpfun(
function(eff,y0.a0,y1.a0,y0.a1,y1.a1,p00,p10,p01,p11,n,m.a0,m.a1) {
 k = length(n) 
 w.n = double()
 w.a0 = double()
 w.a1 = double()
 for (i in 1:k) {
   w.n = c(w.n,rep(1/n[i],n[i]))
   w.a0 = c(w.a0,rep(m.a0[i]/n[i],n[i]))
   w.a1 = c(w.a1,rep(m.a1[i]/n[i],n[i]))
 }  
 if (eff=='DEa0') {
 y0.a0[is.na(y0.a0)==1] = rbinom(length(y0.a0[is.na(y0.a0)==1]),1,p00)
 y1.a0[is.na(y1.a0)==1] = rbinom(length(y1.a0[is.na(y1.a0)==1]),1,p10)
 effect = (1/k)*sum(w.n*(y0.a0-y1.a0))
 }
 if (eff=='DEa1') {
 y0.a1[is.na(y0.a1)==1] = rbinom(length(y0.a1[is.na(y0.a1)==1]),1,p01)
 y1.a1[is.na(y1.a1)==1] = rbinom(length(y1.a1[is.na(y1.a1)==1]),1,p11)
 effect = (1/k)*sum(w.n*(y0.a1-y1.a1))
 }
 if (eff=='IE') {
 y0.a0[is.na(y0.a0)==1] = rbinom(length(y0.a0[is.na(y0.a0)==1]),1,p00)
 y0.a1[is.na(y0.a1)==1] = rbinom(length(y0.a1[is.na(y0.a1)==1]),1,p01)
 effect = (1/k)*sum(w.n*(y0.a0-y0.a1))
 }
 if (eff=='TE') {
 y0.a0[is.na(y0.a0)==1] = rbinom(length(y0.a0[is.na(y0.a0)==1]),1,p00)
 y1.a1[is.na(y1.a1)==1] = rbinom(length(y1.a1[is.na(y1.a1)==1]),1,p11)
 effect = (1/k)*sum(w.n*(y0.a0-y1.a1))
 }
 if (eff=='OE') {
 y0.a0[is.na(y0.a0)==1] = rbinom(length(y0.a0[is.na(y0.a0)==1]),1,p00)
 y1.a0[is.na(y1.a0)==1] = rbinom(length(y1.a0[is.na(y1.a0)==1]),1,p10)
 y0.a1[is.na(y0.a1)==1] = rbinom(length(y0.a1[is.na(y0.a1)==1]),1,p01)
 y1.a1[is.na(y1.a1)==1] = rbinom(length(y1.a1[is.na(y1.a1)==1]),1,p11) 
 effect = (1/k)*sum(w.n*(w.a0*y1.a0+(1-w.a0)*y0.a0-w.a1*y1.a1-(1-w.a1)*y0.a1))
 } 
output.all = list(y0.a0=y0.a0,y1.a0=y1.a0,y0.a1=y0.a1,y1.a1=y1.a1,effect=effect)
return(output.all)
}
)

#Randomizing function
rand = cmpfun(
function(n,m) {
Z = rep(0,n)
trt = sample(n,m)
Z[trt] = 1
return(Z)
}
)

#p-value function
pval = cmpfun(
function(eff,est,null,y0.a0,y1.a0,y0.a1,y1.a1,h,n,m.a0,m.a1,C2) {
k = length(n)
dist=double()
  for (i in 1:C2) {
      g = rand(k,h)
      S = double()
      Z = double()
      w.c = double()
      w.t = double()
      w.n = double()
       for (j in 1:k) {
        S = c(S,rep(g[j],n[j]))
        w.n = c(w.n,rep(1/n[j],n[j]))
        if (g[j]==0) {
        Z = c(Z,rand(n[j],m.a0[j]))
        w.c = c(w.c,rep(1/(n[j]-m.a0[j]),n[j]))
        w.t = c(w.t,rep(1/m.a0[j],n[j]))
        }
        if (g[j]==1) {
        Z = c(Z,rand(n[j],m.a1[j]))
        w.c = c(w.c,rep(1/(n[j]-m.a1[j]),n[j]))
        w.t = c(w.t,rep(1/m.a1[j],n[j]))
        }
      }
      if (eff=='DEa0') {
      ts = (1/(k-h))*sum(w.c*(1-S)*(1-Z)*y0.a0)-(1/(k-h))*sum(w.t*(1-S)*Z*y1.a0)
      dist = c(dist,ts)
      }
      if (eff=='DEa1') {
      ts = (1/h)*sum(w.c*S*(1-Z)*y0.a1)-(1/h)*sum(w.t*S*Z*y1.a1)
      dist = c(dist,ts)
      }
      if (eff=='IE') {
      ts = (1/(k-h))*sum(w.c*(1-S)*(1-Z)*y0.a0)-(1/h)*sum(w.c*S*(1-Z)*y0.a1)
      dist = c(dist,ts)
      }
      if (eff=='TE') {
      ts = (1/(k-h))*sum(w.c*(1-S)*(1-Z)*y0.a0)-(1/h)*sum(w.t*S*Z*y1.a1)
      dist = c(dist,ts)
      }
      if (eff=='OE') {
      ts = (1/(k-h))*sum(w.n*(1-S)*((1-Z)*y0.a0+Z*y1.a0))-(1/h)*sum(w.n*S*((1-Z)*y0.a1+Z*y1.a1))
      dist = c(dist,ts)
    }
    }
    p = mean(abs(dist-null)>=abs(est-null))
return(p)
}
)

#Combination function to compute C1
nchoosem = cmpfun(
function(n,m) {
c = choose(n,m)
trt = combinations(n,m)
Z = matrix(NA,c,n)
for (i in 1:c) {
Z[i,trt[i,]] = 1
Z[i,-trt[i,]] = 0
}
return(Z)
}
)

#Taking care of p_zs on boundary
bd = cmpfun(
  function(x) {
  if (x>1) {x=1}
  if (x<0) {x=0}
  return(x)
}
)

#Locally linear interpolation
lsolve = cmpfun(
  function(x1,y1,x2,y2,level) {
  m = (y2-y1)/(x2-x1)
  return((1/m)*(level-y1)+x1)
}
)

#Compute confidence interval
exactCI = cmpfun(
function(eff,g,data,m.a0,m.a1,B2,C2,level) {
k = length(g)
h = sum(g)
n = apply(data,3,sum)
m = apply(data,3,function(x) x[1,1]+x[1,2])
a = estbound(g,data,m.a0,m.a1)
N = sum(n)

sz00 = sum((1-g)*(n-m))
sz10 = sum((1-g)*m)  
sz01 = sum(g*(n-m))
sz11 = sum(g*m)
  
#compute B1
B1 = double()
if (eff=='DEa0') {B1 = (2^(sz00+sz10))*(4^(N-sz00-sz10))}
if (eff=='DEa1') {B1 = (2^(sz01+sz11))*(4^(N-sz01-sz11))}
if (eff=='IE') {B1 = (2^(sz00+sz01))*(4^(N-sz00-sz01))}
if (eff=='TE') {B1 = (2^(sz00+sz11))*(4^(N-sz00-sz11))}
if (eff=='OE') {B1 = 8^N}

#compute C1
j0 = double()
j1 = double()
for (i in 1:length(n)) {
  j0 = c(j0,choose(n[i],m.a0[i]))
  j1 = c(j1,choose(n[i],m.a1[i]))
}
s = nchoosem(k,h)
nexp = (1-s)*j0+s*j1
C1 = sum(apply(nexp,1,prod))

#targeted search
nhyp = B2/2
prob1 = 1-1/B2
prob2 = 1/B2
count.NA = 0
effect = double()
p = double()
est = double()

if (eff=='DEa0') {
est = a$DEa0.hat
effect = c(effect,a$DEa0.hat)
p = c(p,1)

effect = c(effect,a$DEa0.lb)
p = c(p,pval(eff='DEa0',est=a$DEa0.hat,null=a$DEa0.lb,y0.a0=a$y0.a0.low,y1.a0=a$y1.a0.high,y0.a1=a$y0.a1,y1.a1=a$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

effect = c(effect,a$DEa0.ub)
p = c(p,pval(eff='DEa0',est=a$DEa0.hat,null=a$DEa0.ub,y0.a0=a$y0.a0.high,y1.a0=a$y1.a0.low,y0.a1=a$y0.a1,y1.a1=a$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

l = a$DEa0.hat
lower = a$DEa0.lb  
if(p[effect==a$DEa0.lb]<level) { 
  for (i in 1:nhyp) {
   tl.mid = quantile(c(a$DEa0.lb,l),probs=prob1)
   j = sample.n(eff='DEa0',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.lb+a$y1.a0.ub+tl.mid)/2),p10=bd((a$y0.a0.lb+a$y1.a0.ub-tl.mid)/2),p01=a$y0.a1.hat,p11=a$y1.a1.hat,n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect<=l) {p.test=pval(eff='DEa0',est=a$DEa0.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob1=prob1-1/B2}
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   l = min(effect[p>=level])
 } 
} 

u = a$DEa0.hat
upper = a$DEa0.ub  
if(p[effect==a$DEa0.ub]<level) {
  for (i in 1:nhyp) {
   tu.mid = quantile(c(u,a$DEa0.ub),probs=prob2)
   j = sample.n(eff='DEa0',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.ub+a$y1.a0.lb+tu.mid)/2),p10=bd((a$y0.a0.ub+a$y1.a0.lb-tu.mid)/2),p01=a$y0.a1.hat,p11=a$y1.a1.hat,n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect>=u) {p.test=pval(eff='DEa0',est=a$DEa0.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob2=prob2+1/B2}
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   u = max(effect[p>=level])
 }
}
}

if (eff=='DEa1') {
est = a$DEa1.hat
effect = c(effect,a$DEa1.hat)
p = c(p,1)

effect = c(effect,a$DEa1.lb)
p = c(p,pval(eff='DEa1',est=a$DEa1.hat,null=a$DEa1.lb,y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1.low,y1.a1=a$y1.a1.high,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

effect = c(effect,a$DEa1.ub)
p = c(p,pval(eff='DEa1',est=a$DEa1.hat,null=a$DEa1.ub,y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1.high,y1.a1=a$y1.a1.low,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

l = a$DEa1.hat
lower = a$DEa1.lb  
if(p[effect==a$DEa1.lb]<level) { 
  for (i in 1:nhyp) {
   tl.mid = quantile(c(a$DEa1.lb,l),probs=prob1)
   j = sample.n(eff='DEa1',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=a$y0.a0.hat,p10=a$y1.a0.hat,p01=bd((a$y0.a1.lb+a$y1.a1.ub+tl.mid)/2),p11=bd((a$y0.a1.lb+a$y1.a1.ub-tl.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect<=l) {p.test=pval(eff='DEa1',est=a$DEa1.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob1=prob1-1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   l = min(effect[p>=level])
 }
}

u = a$DEa1.hat
upper = a$DEa1.ub  
if(p[effect==a$DEa1.ub]<level) {
  for (i in 1:nhyp) {
   tu.mid = quantile(c(u,a$DEa1.ub),probs=prob2)
   j = sample.n(eff='DEa1',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=a$y0.a0.hat,p10=a$y1.a0.hat,p01=bd((a$y0.a1.ub+a$y1.a1.lb+tu.mid)/2),p11=bd((a$y0.a1.ub+a$y1.a1.lb-tu.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect>=u) {p.test=pval(eff='DEa1',est=a$DEa1.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob2=prob2+1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]
   u = max(effect[p>=level])
 }
}
}

if (eff=='IE') {
est = a$IE.hat
effect = c(effect,a$IE.hat)
p = c(p,1)

effect = c(effect,a$IE.lb)
p = c(p,pval(eff='IE',est=a$IE.hat,null=a$IE.lb,y0.a0=a$y0.a0.low,y1.a0=a$y1.a0,y0.a1=a$y0.a1.high,y1.a1=a$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

effect = c(effect,a$IE.ub)
p = c(p,pval(eff='IE',est=a$IE.hat,null=a$IE.ub,y0.a0=a$y0.a0.high,y1.a0=a$y1.a0,y0.a1=a$y0.a1.low,y1.a1=a$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

l = a$IE.hat
lower = a$IE.lb  
if(p[effect==a$IE.lb]<level) { 
  for (i in 1:nhyp) {
   tl.mid = quantile(c(a$IE.lb,l),probs=prob1)
   j = sample.n(eff='IE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.lb+a$y0.a1.ub+tl.mid)/2),p10=a$y1.a0.hat,p01=bd((a$y0.a0.lb+a$y0.a1.ub-tl.mid)/2),p11=a$y1.a1.hat,n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect<=l) {p.test=pval(eff='IE',est=a$IE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob1=prob1-1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   l = min(effect[p>=level])
 }
}

u = a$IE.hat
upper = a$IE.ub  
if(p[effect==a$IE.ub]<level) {
  for (i in 1:nhyp) {
   tu.mid = quantile(c(u,a$IE.ub),probs=prob2)
   j = sample.n(eff='IE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.ub+a$y0.a1.lb+tu.mid)/2),p10=a$y1.a0.hat,p01=bd((a$y0.a0.ub+a$y0.a1.lb-tu.mid)/2),p11=a$y1.a1.hat,n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect>=u) {p.test=pval(eff='IE',est=a$IE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}   
   if (is.na(p.test)==1) {prob2=prob2+1/B2}
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   u = max(effect[p>=level])
 }  
}
}

if (eff=='TE') {
est = a$TE.hat
effect = c(effect,a$TE.hat)
p = c(p,1)

effect = c(effect,a$TE.lb)
p = c(p,pval(eff='TE',est=a$TE.hat,null=a$TE.lb,y0.a0=a$y0.a0.low,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1.high,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

effect = c(effect,a$TE.ub)
p = c(p,pval(eff='TE',est=a$TE.hat,null=a$TE.ub,y0.a0=a$y0.a0.high,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1.low,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

l = a$TE.hat
lower = a$TE.lb  
if(p[effect==a$TE.lb]<level) { 
  for (i in 1:nhyp) {
   tl.mid = quantile(c(a$TE.lb,l),probs=prob1)
   j = sample.n(eff='TE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.lb+a$y1.a1.ub+tl.mid)/2),p10=a$y1.a0.hat,p01=a$y0.a1.hat,p11=bd((a$y0.a0.lb+a$y1.a1.ub-tl.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect<=l) {p.test=pval(eff='TE',est=a$TE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob1=prob1-1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   l = min(effect[p>=level])
 }
}

u = a$TE.hat
upper = a$TE.ub  
if(p[effect==a$TE.ub]<level) {
  for (i in 1:nhyp) {
   tu.mid = quantile(c(u,a$TE.ub),probs=prob2)
   j = sample.n(eff='TE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y0.a0.ub+a$y1.a1.lb+tu.mid)/2),p10=a$y1.a0.hat,p01=a$y0.a1.hat,p11=bd((a$y0.a0.ub+a$y1.a1.lb-tu.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect>=u) {p.test=pval(eff='TE',est=a$TE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}   
   if (is.na(p.test)==1) {prob2=prob2+1/B2}
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]   
   u = max(effect[p>=level])
 }
}
}

if (eff=='OE') {
est = a$OE.hat
effect = c(effect,a$OE.hat)
p = c(p,1)

effect = c(effect,a$OE.lb)
p = c(p,pval(eff='OE',est=a$OE.hat,null=a$OE.lb,y0.a0=a$y0.a0.low,y1.a0=a$y1.a0.low,y0.a1=a$y0.a1.high,y1.a1=a$y1.a1.high,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

effect = c(effect,a$OE.ub)
p = c(p,pval(eff='OE',est=a$OE.hat,null=a$OE.ub,y0.a0=a$y0.a0.high,y1.a0=a$y1.a0.high,y0.a1=a$y0.a1.low,y1.a1=a$y1.a1.low,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2))

l = a$OE.hat
lower = a$OE.lb  
if(p[effect==a$OE.lb]<level) { 
  for (i in 1:nhyp) {
   tl.mid = quantile(c(a$OE.lb,l),probs=prob1) 
   j = sample.n(eff='OE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y.a0.lb+a$y.a1.ub+tl.mid)/2),p10=bd((a$y.a0.lb+a$y.a1.ub+tl.mid)/2),p01=bd((a$y.a0.lb+a$y.a1.ub-tl.mid)/2),p11=bd((a$y.a0.lb+a$y.a1.ub-tl.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect<=l) {p.test=pval(eff='OE',est=a$OE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob1=prob1-1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]
   l = min(effect[p>=level])
 }
}

u = a$OE.hat
upper = a$OE.ub  
if(p[effect==a$OE.ub]<level) {
  for (i in 1:nhyp) {
   tu.mid = quantile(c(u,a$OE.ub),probs=prob2)
   j = sample.n(eff='OE',y0.a0=a$y0.a0,y1.a0=a$y1.a0,y0.a1=a$y0.a1,y1.a1=a$y1.a1,p00=bd((a$y.a0.ub+a$y.a1.lb+tu.mid)/2),p10=bd((a$y.a0.ub+a$y.a1.lb+tu.mid)/2),p01=bd((a$y.a0.ub+a$y.a1.lb-tu.mid)/2),p11=bd((a$y.a0.ub+a$y.a1.lb-tu.mid)/2),n=n,m.a0=m.a0,m.a1=m.a1)
   effect = c(effect,j$effect)
   p.test = NA
   if (j$effect>=u) {p.test=pval(eff='OE',est=a$OE.hat,null=j$effect,y0.a0=j$y0.a0,y1.a0=j$y1.a0,y0.a1=j$y0.a1,y1.a1=j$y1.a1,h=h,n=n,m.a0=m.a0,m.a1=m.a1,C2=C2)}
   if (is.na(p.test)==1) {count.NA=count.NA+1}
   if (is.na(p.test)==1) {prob2=prob2+1/B2}   
   p = c(p,p.test)
   effect = effect[is.na(p)==0]
   p = p[is.na(p)==0]
   u = max(effect[p>=level])
 }  
}
}

#Remove NAs; first count fraction
#frac = sum(is.na(p))/length(p)
p = p[is.na(p)==0]
effect = effect[is.na(p)==0]

#Compute intervals
lower1 = min(effect[p>=level])
lower2 = max(effect[p<level & effect<lower1])
upper1 = max(effect[p>=level])
upper2 = min(effect[p<level & effect>upper1])

#Handle missings and NAs if appropriate
if (lower2==-Inf) {
  if (eff=='DEa0') {lower2=a$DEa0.lb}
  if (eff=='DEa1') {lower2=a$DEa1.lb}
  if (eff=='IE') {lower2=a$IE.lb}
  if (eff=='TE') {lower2=a$TE.lb}
  if (eff=='OE') {lower2=a$OE.lb}
}
if (upper2==Inf) {
  if (eff=='DEa0') {upper2=a$DEa0.ub}
  if (eff=='DEa1') {upper2=a$DEa1.ub}
  if (eff=='IE') {upper2=a$IE.ub}
  if (eff=='TE') {upper2=a$TE.ub}
  if (eff=='OE') {upper2=a$OE.ub}
}

#Locally linear interpolation
lower3 = lsolve(lower2,max(p[effect==lower2]),lower1,max(p[effect==lower1]),level)
upper3 = lsolve(upper1,max(p[effect==upper1]),upper2,max(p[effect==upper2]),level)

if (is.na(lower3)==1) {
  if (eff=='DEa0') {lower3=a$DEa0.lb}
  if (eff=='DEa1') {lower3=a$DEa1.lb}
  if (eff=='IE') {lower3=a$IE.lb}
  if (eff=='TE') {lower3=a$TE.lb}
  if (eff=='OE') {lower3=a$OE.lb}
}
if (is.na(upper3)==1) {
  if (eff=='DEa0') {upper3=a$DEa0.ub}
  if (eff=='DEa1') {upper3=a$DEa1.ub}
  if (eff=='IE') {upper3=a$IE.ub}
  if (eff=='TE') {upper3=a$TE.ub}
  if (eff=='OE') {upper3=a$OE.ub}
}

output.all = list(B1=B1,C1=C1,frac.NA=count.NA/B2,prob1=prob1,prob2=prob2,effect=effect,p=p,est=est,lower=lower3,upper=upper3)
return(output.all)
}
)

#k = exactCI('OE',assign.ex,data.ex,rep(3,10),rep(6,10),100,100,0.05)
#plot(k$effect,k$p)

#package.skeleton(name="interferenceCI",code_files="/home/joe/Dropbox/Dissertation/interferenceCI.R")
