gendata.mix <- function(n1,n2,m1,m2,p,alpha,prob.y=c(0.9,0.1))
{
	theta.c1=numeric(p)
	theta.c2=numeric(p)
	
	P=runif(p,0,1)
	
	theta.c1[1]=prob.y[1]
	theta.c2[1]=prob.y[2]
	for(i in 2:p)
        {	theta.c1[i]=rbeta(1,alpha*P[i],alpha*(1-P[i]))
		theta.c2[i]=rbeta(1,alpha*P[i],alpha*(1-P[i]))
	}		
	train=matrix(runif(n1*p+n2*p),n1+n2,p)
	test=matrix(runif(m1*p+m2*p),m1+m2,p)
			
	train=matrix(as.numeric(train<=rbind(matrix(rep(theta.c1,n1),n1,p,TRUE),
            		matrix(rep(theta.c2,n2),n2,p,TRUE))),n1+n2,p)
	test=matrix(as.numeric(test<=rbind(matrix(rep(theta.c1,m1),m1,p,TRUE),
            		matrix(rep(theta.c2,m2),m2,p,TRUE))),m1+m2,p)
	list(train=train,test=test)
}

###############################################################################
indexconvert=function(a)
{
	b=c()
	for(i in 1:length(a))
	{
		if(a[i]==TRUE)          
		b=c(b,i) 
	}
	b
}

