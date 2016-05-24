StabC <-
function(s,p,maxiter=10000,method="monte-carlo")
{
C=0
if(method=="monte-carlo") 
{
  for(i in 1:maxiter) 
  {
    sum_N=0
    N=sample.int(n=s,size=s,replace=TRUE,prob=p)
    for(j in 1:s) 
    {
      sum_N=sum_N+N[j]
      if(sum_N>s) break 
    }
    j=j-1
    C=C+1/j
  }
  C=C/maxiter
}
if(method=="less-memory") 
{
  n=vector(mode="numeric",length=s) 
  for(q in 1:s) 
  {
    parts_q=parts(q)
    num_parts=dim(parts_q)[2]
    for (i in 1:num_parts) 
    {
      part=parts_q[,i]
      for(j in 1:s) n[j]=sum(part==j)
      C=C+factorial(sum(part>0)-1)*prod(p[part])*
        sum(p[(s+1-q):s])/prod(factorial(n))
    }
  }
}
if(method=="more-memory") 
{
  for(q in 1:s) 
  {
    n=array(colSums(array(parts(q),dim=c(q,q*P(q)))==array(rep(1:q,each=q*P(q)),
      dim=c(q,q*P(q)))),dim=c(P(q),q))
    probs=t(array(p[1:q],dim=c(q,P(q))))^n/factorial(n)
    C=C+sum(factorial(rowSums(n)-1)*apply(probs,1,prod)*sum(p[(s+1-q):s]))
  }
}
if(method=="progressive") 
{
  n=vector(mode="numeric",length=s) 
  for(q in 1:s) 
  {
    parts_q=parts(q)
    num_parts=dim(parts_q)[2]
    for (i in 1:num_parts) 
    {
      part=parts_q[,i]
      #for(j in 1:s) n[j]=sum(part==j)
      C=C+multinom(part[part>0])*prod(p[part])*sum(p[(s+1-q):s])/sum(part>0)
    }
  }
}
return(1/C)
}
