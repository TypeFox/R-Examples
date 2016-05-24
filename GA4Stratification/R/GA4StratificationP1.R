GA4StratificationP1 <-
function(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)
{
   dataName=data.frame(dataName)
   dataName=data.frame(dataName[order(dataName[,1]),])
   tableData=as.data.frame(table(dataName))
   randomnumRange=cumsum(tableData[,2])
   lengthRandomnum=length(randomnumRange)

   lengthData=nrow(dataName)
   randomGeneration=array(0,dim=c(GAgenerationSize,lengthData))

   nocrom=nrow(randomGeneration)
   fitp1=array(0,dim=c(1,nocrom))
   fit=array(0,dim=c(nocrom,1))
   N=means=s=n=vars=mas=NN=k=p=t=array(0,dim=c(nocrom,numberOfStrata))
   dd=array(0,dim=c(nocrom,1))

   randomNumbers=array(0,dim=c(GAgenerationSize,numberOfStrata-1))

   randomNumbersX=array(0,dim=c(nocrom,2))

   aftermut=array(0,dim=c(GAgenerationSize,lengthData+numberOfStrata+1))
   cumTotal=cumsum(dataName)
   sumSquares=cumsum(dataName^2)


   for (i in 1:GAgenerationSize)
   {
	randomNumbers[i,]=randomnumGenerator(randomnumRange,lengthRandomnum,numberOfStrata-1)

   }
   

   son=array(lengthData,dim=c(GAgenerationSize,1))
   indis=array(c(1:GAgenerationSize,randomNumbers,son),dim=c(GAgenerationSize,(numberOfStrata+1)))

   for(i in 2:(numberOfStrata+1))
   {
   	randomGeneration[indis[,c(1,i)]]=1
   }

   bestValue=-99999999999999999999999999999999999999999999999999999999999999
   for ( i in 1:iteration )
   {

	fitnessValueGeneration=GA4StratificationP1fit(randomGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
	if ( max(fitnessValueGeneration)>bestValue )
	{
	    bestValue=max(fitnessValueGeneration)
	    bestGeneration=randomGeneration[max(which(fitnessValueGeneration==bestValue)),]
	}

	randomGeneration=GA4StratificationSelection(randomGeneration,fitnessValueGeneration)

	randomGeneration=GA4StratificationP1x(randomGeneration,bestGeneration,dataName,numberOfStrata,sampleSize,fitnessValueGeneration,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t,randomNumbersX,tableData,randomnumRange,lengthRandomnum)

	randomGeneration=GA4StratificationP1m(randomGeneration,mutationRate,i)

	fitnessValueGeneration=GA4StratificationP1fit(randomGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)

	if ( max(fitnessValueGeneration)>bestValue )
	{
	    bestValue=max(fitnessValueGeneration)
	    bestGeneration=randomGeneration[max(which(fitnessValueGeneration==bestValue)),]
	} else

	{
	randomGeneration[sample(GAgenerationSize,1),]=bestGeneration
	}

	cat(i, " ",-bestValue,'\n')
	flush.console()
	Sys.sleep(1)
   }
   GA4StratificationP1fitt(bestGeneration,dataName,numberOfStrata,sampleSize,-bestValue,cumTotal,sumSquares)
}

