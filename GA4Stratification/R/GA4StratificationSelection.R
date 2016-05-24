GA4StratificationSelection <-
function(selectionGeneration,selectionGenerationFitness)
   {
   	rowSelectionGeneration=nrow(selectionGeneration)
	colSelectionGeneration=ncol(selectionGeneration)
   	selectionStrata=array(0,dim=c(rowSelectionGeneration,(colSelectionGeneration+1)))
   	newSelectionGeneration=cbind(selectionGeneration, -selectionGenerationFitness)
   	sortedSelectionGeneration=newSelectionGeneration[order(newSelectionGeneration[,colSelectionGeneration+1]),]
   	sortedSelectionGenerationFitness=selectionGenerationFitness[order(-selectionGenerationFitness),]
   	wheelOld = sortedSelectionGenerationFitness / sum(selectionGenerationFitness)
	wheel=1:rowSelectionGeneration

	for(i in rowSelectionGeneration:1)
	{
		wheel[rowSelectionGeneration+1-i]=wheelOld[i]
	}
	wheel=cumsum(wheel)

   	for ( i in  1:rowSelectionGeneration )
   	{
		r = runif(1,0,1)
		for ( j in  1:rowSelectionGeneration )
		{
			if(r < wheel[j])
			{
				selectionStrata[i,] = sortedSelectionGeneration[j,]
				break;
			}
		}
	}
	randomGeneration=selectionStrata[,1:colSelectionGeneration]
	fitnessValueGeneration=selectionStrata[,(colSelectionGeneration+1)]
	return(randomGeneration)
   }

