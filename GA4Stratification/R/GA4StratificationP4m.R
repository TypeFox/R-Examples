GA4StratificationP4m <-
function(mutationGeneration,numberOfStrata,mutationRate,it)
   {

   	rowOfMutationGeneration=nrow(mutationGeneration)
   	colOfMutationGeneration=ncol(mutationGeneration)
	c=colOfMutationGeneration-numberOfStrata
   	for ( i in 1:rowOfMutationGeneration)
   	{
    	   if ( runif(1,0,1) < mutationRate )
 	   {
		for ( k in 1:2 )
  		{
	         if ( runif(1,0,1) < 0.35 )
		   {
			if ( it<100 )
			{
	               ones=which(mutationGeneration[i,1:c]==1,arr.ind=TRUE)
      	         zeros=which(mutationGeneration[i,1:c]==0,arr.ind=TRUE)
            	   mutationPoint=ones[sample(1:(length(ones)-1),1)]
	               mutationPoint1=zeros[sample(1:length(zeros),1)]
	               mutationGeneration[i,mutationPoint]=0
	               mutationGeneration[i,mutationPoint1]=1
            	} else
			{ 
	               ones=which(mutationGeneration[i,1:c]==1,arr.ind=TRUE)
            	   mutationPoint=ones[sample(1:(length(ones)-1),1)]
                     if ( runif(1,0,1)<0.51 )
		         {
				if ( mutationGeneration[i,(mutationPoint+1)]==0 )
				{
					mutationGeneration[i,mutationPoint]=0
		    			mutationGeneration[i,(mutationPoint+1)]=1
              		}
			   } else if ( mutationPoint>1)
			   {
				if ( mutationGeneration[i,(mutationPoint-1)]==0 ) 
				{
				      mutationGeneration[i,mutationPoint]=0
				      mutationGeneration[i,mutationPoint-1]=1    
				}
			   }
			}
    		} else
		{
		   mutationPoint=sample((c+1):(c+numberOfStrata),1)
		   mutationPoint1=sample((c+1):(c+numberOfStrata),1)

	         while  ( mutationPoint == mutationPoint1 )
		   {
			mutationPoint=sample((c+1):(c+numberOfStrata),1)
		   	mutationPoint1=sample((c+1):(c+numberOfStrata),1)
		   }
     
               if (mutationGeneration[i,mutationPoint]>2 & mutationGeneration[i,mutationPoint1]>2)
		   {
			if (runif(1,0,1)< 0.51)
			{
				mutationGeneration[i,mutationPoint]=mutationGeneration[i,mutationPoint]-1
		      	mutationGeneration[i,mutationPoint1]=mutationGeneration[i,mutationPoint1]+1

           		} else
			{
		            mutationGeneration[i,mutationPoint]=mutationGeneration[i,mutationPoint]+1
		            mutationGeneration[i,mutationPoint1]=mutationGeneration[i,mutationPoint1]-1           
			}

 		   } else if (mutationGeneration[i,mutationPoint]==2 & mutationGeneration[i,mutationPoint1]>2)
		   {
			mutationGeneration[i,mutationPoint]=mutationGeneration[i,mutationPoint]+1
           		mutationGeneration[i,mutationPoint1]=mutationGeneration[i,mutationPoint1]-1

          	   } else (mutationGeneration[i,mutationPoint]>2 & mutationGeneration[i,mutationPoint1]==2)

	            mutationGeneration[i,mutationPoint]=mutationGeneration[i,mutationPoint]-1
       	      mutationGeneration[i,mutationPoint1]=mutationGeneration[i,mutationPoint1]+1
		   }
		}
	  }
	}
   return(mutationGeneration)
  }

