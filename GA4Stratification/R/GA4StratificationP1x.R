GA4StratificationP1x <-
function(crossoverGeneration,bestGeneration,dataName,numberOfStrata,sampleSize,fitnessValueGeneration,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t,randomNumbersX,tableData,randomnumRange,lengthRandomnum)
 {

   fitnessValueParents=fitnessValueGeneration
   parents=cbind(crossoverGeneration,fitnessValueParents)
   crossoverGenerationp=crossoverGeneration
   rowCrossoverGenerationp=nocrom
	   
   for (i in 1:rowCrossoverGenerationp)
   {
	randomNumbersX[i,]=randomnumGenerator((1:rowCrossoverGenerationp),(rowCrossoverGenerationp+1),2)
   }

   mother=father=NULL

   for (i in 1:rowCrossoverGenerationp)
   {
	  mother=randomNumbersX[i,1]
	  father=randomNumbersX[i,2]
   
    	  crossoverPoint=sample((randomnumRange[1:lengthRandomnum-1]),1)
        while ( sum(crossoverGenerationp[mother,c(1:crossoverPoint)])!=sum(crossoverGenerationp[father,c(1:crossoverPoint)]) )
        {
		crossoverPoint=sample((randomnumRange[1:lengthRandomnum-1]),1)

	  }
    	  crossoverGeneration[i,c(1:crossoverPoint)]=crossoverGenerationp[mother,c(1:crossoverPoint)]
        crossoverGeneration[i,c((crossoverPoint+1):lengthData)]=crossoverGenerationp[father,c((crossoverPoint+1):lengthData)]
   }
      s=GA4StratificationP1fit(crossoverGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
      crossoverGenerationx=cbind(crossoverGeneration,s)
      GA4StratificationP1x=rbind(parents, crossoverGenerationx)
      GA4StratificationP1x=GA4StratificationP1x[order(GA4StratificationP1x[,(lengthData+1)]),]
      GA4StratificationP1x=GA4StratificationP1x[c((rowCrossoverGenerationp+1):(rowCrossoverGenerationp*2)),c(1:lengthData)]
	return(GA4StratificationP1x)
 }

