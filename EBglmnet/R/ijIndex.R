ijIndex <-function(trueLoc,K){
#for easier comparison, we write a function to convert the epistatic effect into interaction pairs:
	nEff = length(trueLoc);
	index = matrix(0,0,2);
	Used1 = trueLoc;
	for (i in 1:K)
	{ 						
		Used1 = Used1-(K-i+1);
		I2 = which(Used1<=0);
		if(length(I2)>0)
		{
			Used2 = Used1[I2];                            
			temp = length(Used2);                        
			if (i == 1)
			{
				locus1  = Used2 + K-i+1;                        
				locus2  = locus1;
			}else
			{
				locus1 = (i -1)*rep(1,temp);                  
				locus2  = Used2 +K;                             
			}
			index  = rbind(index,cbind(locus1, locus2));
			Used1  = Used1[-I2];
		}	
	}	
	return(index);
}#end of function definition.

