var.e<-
function(obs.t,delta,z,r,dx,iter.max,n,p,num.sim){
  beta.est = matrix(nrow=num.sim,ncol=p)
  Rt.est <-Rt0.est <- matrix(nrow=num.sim,ncol=n)
  for(i in 1:num.sim){
	wt = rexp(n)
	beta.ini = rep(0,p)
	Rt.ini = cumsum(1/(n:1))
	temp = solve.beta(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max)
	beta.est[i,] = temp$Beta
	Rt.est[i,] = temp$data$Rt
      if(p>1) Rt0.est[i,]<-temp$data$Rt*exp(-sum(temp$Beta*apply(z,2,mean)))
      if(p==1) Rt0.est[i,]<-temp$data$Rt*exp(-sum(temp$Beta*mean(z)))
  }
return(list(beta.est,Rt0.est))
}
