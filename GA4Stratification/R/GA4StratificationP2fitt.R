GA4StratificationP2fitt <-
function(bestGeneration,dataName,numberOfStrata,sampleSize,bestValue,cumTotal,sumSquares)
{
   c=nrow(dataName)
   nocrom=length(bestGeneration)/c
   bestGeneration=array(bestGeneration,dim=c(length(bestGeneration)/c,c))
   fitp1=array(0,dim=c(1,nocrom))
   fit=array(0,dim=c(nocrom,1))
   N=means=s=n=vars=mas=NN=k=p=t=array(0,dim=c(nocrom,numberOfStrata))

   dd=array(0,dim=c(nocrom,1))

   for ( i in 1:nocrom )
   {

   	mas[i,]=which(bestGeneration[i,]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5

    	for ( j in 2:numberOfStrata )
		{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
			N
	      	means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
			means
      		s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
			s
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
        	p[i,]=which(k[i,]==t[i,],arr.ind=TRUE)
        	n[i,p[i,1]]=min(n[i,p[i,1]]+sampleSize-sum(n[i,]),N[i,p[i,1]])
		}
   }

   return(array(c(N,n,bestValue),dim=c(numberOfStrata,3)))
}

