"getJac" <- function (k,x) 
{
  jac<-array(0,dim=c(length(x),length(k),length(k)))

  for(i in 1:length(k))
	jac[,i,i]<- -x * exp(-k[i]*x)

  jac

}