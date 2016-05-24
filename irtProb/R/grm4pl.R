`grm4pl` <-
function(N=10,theta=0,S=0,C=0,D=0,s=1/1.702,b=0,c=0,d=1) {
 rbinom(N,1,pm4pl(theta=theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d))
 }
