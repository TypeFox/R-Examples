DTMC <-
function(tmat,io,N,trace) 
{
sq<-ifelse(nrow(tmat)==ncol(tmat),1,0)
sum<-ifelse(sum(tmat)==nrow(tmat),1,0)
pp<-ifelse(sum(io)==1 & length(io)==nrow(tmat),1,0)
stopifnot(sq+sum+pp==3)
nr<-nrow(tmat)  
p<-round(upper.tri(matrix(1, nrow(tmat), nrow(tmat)), diag = TRUE))
cs<-tmat%*%p   
U<-runif(N,0,1) 
U1<-runif(1,0,1)
cp<-c(0,cumsum(io)) # I use this formulation to create a cumulative sum to condense the number of logical operations required to compute each iteration.
s0=t(ifelse(cp[2:(nr+1)]>=U1,1,0)*ifelse(cp[1:nr]<U1,1,0)) # this is to determine the initial state
s=matrix(s0,ncol=1) 

state=matrix(ncol=N,nrow=nrow(tmat))

for(i in 1:N)
{
  state[,i]=s # first iteration of S defined above
  ppi=c(0,t(s)%*%cs)
  s=matrix(t(ifelse(ppi[2:(nr+1)]>U[i],1,0)*ifelse(ppi[1:nr]<=U[i],1,0)),ncol=1)
  
}
# We can determine long run probabilities of each state through the WLLN, but this method is not generally computationally tractable.
# Therefore, it is commented out
#print(sum(state[1,])/N)
#print(sum(state[2,])/N)
#print(sum(state[3,])/N)
#print(state)
if(trace==TRUE)
{
  
  gg<-apply(state,2,which.max)
  
 print( plot(gg,type="l",ylab="State",xlab="Number of Iterations",yaxp=c(1,nrow(tmat),nrow(tmat)-1)))
}
return(state)
}

