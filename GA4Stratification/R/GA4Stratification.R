GA4Stratification <-
function(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate,sampleAllocation)
{
   if(sampleAllocation == "Equal")
   {
	GA4StratificationP1(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)
	
   } else if(sampleAllocation=="Proportional") {
	
	GA4StratificationP2(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)

   } else if(sampleAllocation == "Neyman") {
	
	GA4StratificationP3(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)
   
   } else if (sampleAllocation == "GA") {

	GA4StratificationP4(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)

   }

}

