#####Defining initial mu and var
Initial_mu_var=function(rstat,wholeindex)
{
  upmu<-rep(0,2)
  upvar<-rep(0,2)
  initial_mu_var=list()
  initial_mu_var$upmu[1]=mean(rstat[which(wholeindex==0)])
  initial_mu_var$upmu[2]=mean(rstat[which(wholeindex==1)])
  initial_mu_var$upvar[1]=var(rstat[which(wholeindex==0)])
  initial_mu_var$upvar[2]=var(rstat[which(wholeindex==1)])
  return(initial_mu_var)
}
