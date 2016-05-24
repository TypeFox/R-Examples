X.pathway=function(n,p,q1,q2){
if(q1+q2>p){warning("q1+q2>p")}
#q1=sum(beta>0) ###  the number of beta>0 ###
#q2=sum(beta<0) ###  the number of beta<0 ###
q0=p-q1-q2 ### the number of beta=0 ###

X=matrix(0,n,p)
for(i in 1:n){
  if(q0>0){
    X[i,(q1+q2+1):p]=runif(q0,min=-1.5,max=1.5)/0.866
  }
  if(q1>0){
    A1=runif(1,min=-0.75,max=0.75)
    X[i,1:q1]=( runif(q1,min=-0.75,max=0.75)+A1 )/0.612
  }   
  if(q2>0){
    A2=runif(1,min=-0.75,max=0.75)
    X[i,(q1+1):(q1+q2)]=( runif(q2,min=-0.75,max=0.75)+A2 )/0.612
  }
}
X
}
