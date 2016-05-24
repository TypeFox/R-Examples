GA4StratificationP4fitt <-
function(bestGeneration,dataName,numberOfStrata,sampleSize,bestValue,cumTotal,sumSquares)
{
   c=nrow(dataName)
   nocrom=length(bestGeneration)/(c+numberOfStrata)
   bestGeneration=array(bestGeneration,dim=c(length(bestGeneration)/(c+numberOfStrata),(c+numberOfStrata)))
   fitp1=array(0,dim=c(1,nocrom))
   fit=array(0,dim=c(nocrom,1))
   N=means=s=n=vars=mas=NN=k=p=t=array(0,dim=c(nocrom,numberOfStrata))
   dd=array(0,dim=c(nocrom,1))

   for ( i in 1:nocrom )
   {

   	mas[i,]=which(bestGeneration[i,1:c]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5
      n[i,]=bestGeneration[i,(c+1):(c+numberOfStrata)]

    	for ( j in 2:numberOfStrata )
	{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
			N
	      	means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
			means
      		s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
			s
	}   

   }

   return(array(c(N,n,bestValue),dim=c(numberOfStrata,3)))
}

