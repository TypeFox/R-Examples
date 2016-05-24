retlevel.gpd.graph <-
function(vector,data,threshold)
  {
  Nu=length(data[data>threshold])
  N=length(data)
  plimiar=1-Nu/N
  t=seq(trunc(1/(1-plimiar)+1),100,1)
  n=length(t)
  li=array(0,c(n))
  ls=array(0,c(n))
  pred=array(0,c(n))
  for (s in 1:n)
        {amostra=qgpd(1-(1/(s+min(t)))*N/Nu,vector[,2],threshold,vector[,1])
         li[s]=quantile(amostra,0.025)   
         ls[s]=quantile(amostra,0.975)
         pred[s]=quantile(amostra,0.5)
        }
  plot(t,pred,type="l",ylim=c(max(min(li),0),max(ls)),ylab="returns")
  lines(t,li,lty=2)
  lines(t,ls,lty=2)
  }
