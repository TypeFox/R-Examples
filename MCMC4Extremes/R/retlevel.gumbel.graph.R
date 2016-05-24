retlevel.gumbel.graph <-
function(vector)
  {
  t=seq(1,100,1)
  n=length(t)
  li=array(0,c(n))
  ls=array(0,c(n))
  pred=array(0,c(n))
  for (s in 1:n)
        {amostra=vector[,1]-vector[,2]*log(-log(1-1/s))
         li[s]=quantile(amostra,0.025)   
         ls[s]=quantile(amostra,0.975)
         pred[s]=quantile(amostra,0.5)
        }
  plot(t,pred,type="l",ylim=c(li[2],max(ls)),ylab="returns")
  lines(t,li,lty=2)
  lines(t,ls,lty=2)
  li
  }
