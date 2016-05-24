GA4StratificationP1m <-
function(mutationGeneration,mutationRate,it)
   {

   	rowOfMutationGeneration=nrow(mutationGeneration)
   	for ( i in 1:rowOfMutationGeneration )
   	{
   	   for ( k in 1:5 )
  	   {
    		if ( runif(1,0,1) < mutationRate )
 			{
				if ( runif(1,0,1) < mutationRate )
				{
					if ( it<50 )
			 		{
	               		ones=which(mutationGeneration[i,]==1,arr.ind=TRUE)
      	         		zeros=which(mutationGeneration[i,]==0,arr.ind=TRUE)
            	   		mutationPoint=ones[sample((length(ones)-1),1)]
	               		mutationPoint1=zeros[sample(length(zeros),1)]
	               		mutationGeneration[i,mutationPoint]=0
	               		mutationGeneration[i,mutationPoint1]=1
            	 	} else
			 		{ 
	               		ones=which(mutationGeneration[i,]==1,arr.ind=TRUE)
            	   		mutationPoint=ones[sample((length(ones)-1),1)]
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
							       mutationGeneration[i,mutationPoint]=0;
							       mutationGeneration[i,mutationPoint-1]=1;     
							}
						 }
					}
    			}
			}
  		}
   	}
   return(mutationGeneration)
   }

