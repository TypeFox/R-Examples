allele.sharing <- function(pop1, pop2, allele.column, onlypairs=FALSE, value=NA, ref.pop)
	{ 

	# Function calculates sharing rate of alleles for pop1 and pop2 for each locus in column [allele.column] and [allele.column+1]
	# Output:  tab.all == matrix 
  	# rows=pop1; col=pop2
	#
	# 	ind1 	ind2 	ind3 	ind4
	# ind1 	4 	0 	2 	3	 
	# ind2 	4 	0 	2 	3
	# ind3 	4 	0 	2 	3
	# .    	. 	. 	. 	.
  
	# Preparing matrix	
	matrix.share <- matrix(numeric(0),nrow=length(pop1[,1]),ncol=length(pop2[,1]))
	row.names(matrix.share) <- pop1[,1]
	colnames(matrix.share) <- pop2[,1]
	pop.size <- length(pop1[,1])+length(pop2[,1])-2

	if (onlypairs==FALSE)
		
		{

			
	for (i in 1:(length(matrix.share)))

		{


			col.position <- ceiling(i/length(row.names(matrix.share)))  #rounds up
			row.position <- (i-(length(row.names(matrix.share))*(ceiling(i/length(row.names(matrix.share)))-1)))
			
			if (col.position<row.position)
			{

	
  	# Estimate compare Li and Horvitz 1953
  	if (value=="Bxy") {
      matrix.share[row.position,col.position]<-sum(
	as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]),
	as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1]),
	as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]),
	as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1]==pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1])
	)/4	}

	# Estimate compare Blouin 1996
	if (value=="Mxy") {
		matrix.share[row.position,col.position] <- 
					
						sum(c(
							as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column]),
							as.numeric(pop1[which(pop1[,1]==row.names(matrix.share)[row.position]),allele.column+1])
							)== 
						c(
							as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column]),
							as.numeric(pop2[which(pop2[,1]==colnames(matrix.share)[col.position]),allele.column+1])
							))/2
					
			}
	
  # Estimate compare Queller and Goodnight 1989 and Oliehoek et al. 2006
	if (value=="rxy") 
		{  
    
   	matrix.share[row.position,col.position] <- queller(pop1,pop2,allele.column,row.position,col.position,matrix.share,ref.pop)
			
	  }
		}
		}


		
}
	if (onlypairs==TRUE)
	{
	     
	     for (i in 1:length(colnames(matrix.share)))
        	{

			# Estimate compare Li and Horvitz 1953
			if (value=="Bxy") 
				{
         matrix.share[i,i] <- sum(mean(as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column][1]==pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column][1])),mean(as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column][1]==pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column+1][1])),mean(as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column+1][1]==pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column][1])),mean(as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column+1][1]==pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column+1][1])))/4
				}

			# Estimate compare Blouin 1996		
			if (value=="Mxy") 
				{
				matrix.share[i,i] <- 
					sum(c(
							as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column]),
							as.numeric(pop1[which(pop1[,1]==row.names(matrix.share[i,i,drop=FALSE])),allele.column+1])
							)
							== 
					    c(	
							as.numeric(pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column]),
							as.numeric(pop2[which(pop2[,1]==colnames(matrix.share[i,i,drop=FALSE])),allele.column+1])
							))/2
				}

			# Estimate compare Queller and Goodnight 1989		
			if (value=="rxy") 
				{
         
				matrix.share[i,i] <- queller(pop1,pop2,allele.column,i,i,matrix.share,ref.pop)        
        
				}
        	}

		
	
		}

return(matrix.share)		

	}
