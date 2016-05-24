offspring <- function(parent1,parent2,allele.column, pairs)
		{
		
			# Offspring is calculated from parent1 crossed with parent2 (class = input/data.frame)
			# Output == population according to format from Input.txt; offspr (class = data.frame)
			# allele.column... allele.column+1 as source for offspring
		
		individual <- vector(mode="character")
		population <- vector(mode="character")
		locusA <- vector(mode="numeric")
		locusB <- vector(mode="numeric")
		offspring <- data.frame(individual,population,locusA,locusB)

		for (k in 1:pairs)
			
			{
      # Mendel
			off1 <- c(parent1[k,allele.column],parent2[k,allele.column])
			off2 <- c(parent1[k,allele.column],parent2[k,allele.column+1])
			off3 <- c(parent1[k,allele.column+1],parent2[k,allele.column])
			off4 <- c(parent1[k,allele.column+1],parent2[k,allele.column+1])
			offall <- data.frame(off1,off2,off3,off4)
			
			
      # Random samples
			
				a <- sample(seq(1:4),1)
				individual <- paste(parent1[k,1],parent2[k,1],"-",k,"-",a,sep="")
				population <- paste("off_",parent1[1,2],parent2[1,2],sep="")
				locusA <- min(offall[,a])
				locusB <- max(offall[,a])
				offspring <- rbind(offspring,(data.frame(individual,population,locusA,locusB)))
				
			}   
			
			offspring <- remove.na.rows(offspring)
			return(offspring)
			
		}
		
