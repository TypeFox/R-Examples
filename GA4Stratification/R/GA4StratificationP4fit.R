GA4StratificationP4fit <-
function(randomGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,c,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
{

   for ( i in 1:nocrom )
   {
   	mas[i,]=which(randomGeneration[i,1:c]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
      if(N[i,1]==1)
	{s[i,1]=0
	} else
	{
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5
	}
      n[i,]=randomGeneration[i,(c+1):(c+numberOfStrata)]

    	for ( j in 2:numberOfStrata )
	{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
		means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
		means
      	if(N[i,j]==1)
		{s[i,j]=0
		} else
		{
      	   s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
		}
	}   

	for ( j in 1:numberOfStrata )
	{
         	vars[i,j]=((N[i,j]-n[i,j])*s[i,j]^2*N[i,j]^2)/(c^2*n[i,j]*N[i,j])
	}

    	dd[i,]=min((N[i,]-n[i,]))
    	NN[i,]=cumsum(N[i,])
	fit[i,]=sum(vars[i,])
	
    	if ( dd[i]<0 ) 
	{
		fit[i]= 99999999999999999999999
    	} else if (!all(N[i,]!=1)) 
	{
		fit[i]= 99999999999999999999999
    	} else if (!all(N[i,]!=0)) 
	{
		fit[i]= 99999999999999999999999
	} else
	{
		fit[i]=fit[i]
	}
   

    	p4fit=array(-fit,dim=c(nocrom,1))

   }
   return(p4fit)
}

