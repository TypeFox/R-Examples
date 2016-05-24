# file OrdinalLogisticBiplot/R/plotOrdinalFittedVariable.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


plotOrdinalFittedVariable <- function(nameVariable,coeffic,D,numFactors,planex = 1,planey = 2,xi=-3.5,xu=3.5,yi=-3.5,yu=3.5,
        margin=0,CexVar=0.7,ColorVar="blue",PchVar=0.7,addToPlot=FALSE,showIIC = TRUE,iicxi=-2.5,iicxu=2.5){

   numcat = length(coeffic) - numFactors + 1
   
   if(D == 1.702){
      coeffic = c(coeffic[1:numFactors],(-1)*coeffic[(numFactors+1):length(coeffic)])
   }   
   
   ordBipVar = CalculateOrdinalBiplotGeneral(nameVariable,numcat,coeffic,planex,planey,numFactors,D)

   if(addToPlot == FALSE){
      dev.new()
      plot(0, 0, cex = 0,asp=1, xaxt = "s", yaxt = "s" ,xlim=c(xi,xu),ylim=c(yi,yu),
          main="Ordinal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))
   }
   plot.ordBipVariable(ordBipVar,D,planex,planey,xi,xu,yi,yu,margin,numFactors,CexVar,ColorVar,PchVar)
   
   if(showIIC){
      coeffic = ordBipVar$coef
      slopeort = ordBipVar$slope
      plotCurvesCategoriesVariable(coeffic,slopeort,D,numcat,nameVariable,iicxi,iicxu,planex,planey)
   }

}
