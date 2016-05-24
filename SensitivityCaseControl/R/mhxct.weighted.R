mhxct.weighted <-
function(mi,Ji,gamma=1,m=1){
#Cases (presumably narrow) getting a weight of m
#Exact Mantel-Haenszel-Birch-Cox test for case-control
#sensitivity analysis
#If Ji is a number, there are Ji subjects in every set
#If Ji is a vector, set i has Ji[i] subjects
#Each set contains one case.  There are mi[i] exposed
#subjects in set i
I<-length(mi)
if (length(Ji)==1) Ji<-rep(Ji,I)
for (i in 1:I){
pi<-gamma*mi[i]/((gamma*mi[i])+(Ji[i]-mi[i]))
if (i==1){
g<-c(1-pi,rep(0,m-1),pi)
}
  else {g<-gconv(g,c(1-pi,rep(0,m-1),pi))}
}
   g
}

