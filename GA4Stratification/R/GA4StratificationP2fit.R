GA4StratificationP2fit <-
function(randomGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,c,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
{
   for ( i in 1:nocrom )
   {
   	mas[i,]=which(randomGeneration[i,]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
      if(N[i,1]==1)
	{s[i,1]=0
	} else
	{
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5
	}
    	for ( j in 2:numberOfStrata )
	{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
      	means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
      	if(N[i,j]==1)
		{s[i,j]=0
		} else
		{

      		s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
		}
	}   

    	for ( j in 1:numberOfStrata )
		{
      		n[i,j]=max(1,floor(sampleSize*N[i,j]/sum(N[i,])))
        	n[i,j]=min(n[i,j],N[i,j])
		}
		if ( sampleSize-sum(n[i,])>0 )
		{
      	k[i,]=N[i,]-n[i,]
        	t[i,]=max(k[i,])
        	p[i,]=which(k[i,]==t[i,],arr.ind=TRUE)[1]
        	n[i,p[i,1]]=min(n[i,p[i,1]]+sampleSize-sum(n[i,]),N[i,p[i,1]])
		}

		for ( j in 1:numberOfStrata )
		{
        	vars[i,j]=((N[i,j]-n[i,j])*s[i,j]^2*N[i,j]^2)/(c^2*n[i,j]*N[i,j])
		}

    	dd[i,]=min((N[i,]-n[i,]))
    	NN[i,]=cumsum(N[i,])
    	kl=0

   		fit[i,]=sum(vars[i,])

    	if ( dd[i]<0 ) 
		{
			fit[i]= 9999999999999999
    	} else if (!all(N[i,]!=1)) 
		{
			fit[i]= 999999999999999999
    	} else if (!all(N[i,]!=0)) 
		{
			fit[i]= 999999999999999999
		} else
		{
			fit[i]=fit[i]
		}
   
    	for ( j in 1:(numberOfStrata-1) )
		{
        	kl=kl+dataName[(NN[i,j]+1),]-dataName[NN[i,j],]
		}

    	p2fit=array(-fit,dim=c(nocrom,1))

   }
   return(p2fit)
}

