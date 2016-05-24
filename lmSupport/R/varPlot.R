varPlot <-
function(TheVar, VarName = '', IDs = NULL, AddPoints = 'Strip', AddDensity=TRUE, Detail=2)
#creates a histogram (using default parameters for hist), with options to include
#rug plot, density plot and univariate descriptives of varying detail
#Revision history
#2011-04-11, fixed bug with double plot of hist, JJC
#2011-08-23:  added ID function with default of no ID, JJC
{
  #print descriptives
  print(paste('Descriptive statistics: ',VarName ,sep=''))
  print(varDescribe(TheVar,Detail))
    
  par(cex.lab=1.5, cex.axis=1.2, lwd=2)  #set graphics parameters for pretty graphs
  HistData = hist(TheVar, main='', xlab=VarName)
  
  switch(toupper(AddPoints),
     STRIP = {figStripChart(TheVar,col='red', , cex= 0.3)},
     RUG = {rug(TheVar,col='red')}
  )

  if(AddDensity) 
  {
      DensityData = density(TheVar)
      lines(DensityData$x,DensityData$y * (max(HistData$counts)/max(DensityData$y)), col='blue')  
  }
  
  if (!is.null(IDs)) 
  {
     Indices = identify(TheVar,rep(0,length(TheVar)), labels=IDs)
     return(Cases = list(Indices=Indices, Rownames = IDs[Indices], Values = TheVar[Indices]))
  } 
  
}