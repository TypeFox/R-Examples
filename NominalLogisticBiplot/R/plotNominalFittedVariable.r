# file NominalLogisticBiplot/R/plotNominalFittedVariable.R
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

plotNominalFittedVariable <- function(nameVar,numcateg,beta,varstudyC,rowCoords,levelsVar=NULL,numFactors=2,planex = 1,planey = 2,xi=-3.5,xu=3.5,yi=-3.5,yu=3.5,
        CexVar=0.7,ColorVar="blue",PchVar=0.7,addToPlot=FALSE,QuitNotPredicted=TRUE,ShowResults=TRUE,
        linesVoronoi=TRUE,LabelVar=TRUE){
   
   mvv = mvvSingleVariable(nameVar,numcateg,beta,varstudyC,rowCoords,planex,planey,numFactors,
                          QuitNotPredicted,ShowResults)
                           
   AtLeastR2 = 0.01
   LabValVar=levelsVar 
   if(addToPlot == FALSE){
      dev.new()  
      plot(0, 0, cex = 0,asp=1, xaxt = "s", yaxt = "s" ,xlim=c(xi,xu),ylim=c(yi,yu),
          main="Nominal Logistic Biplot", xlab=paste("Axis ",planex,sep=""), ylab=paste("Axis ",planey,sep=""))      
   }
   
   if(mvv$numFit==1){
        Barx=sum(rowCoords[,1])/nrow(rowCoords)
        Bary=sum(rowCoords[,2])/nrow(rowCoords)
        points(Barx,Bary,pch=PchVar,cex=CexVar,col=ColorVar)
        
        text(Barx,Bary, paste(mvv$equivFit,"_",nameVar ,sep="") , col = ColorVar, cex = CexVar,pos=1,offset=0.1)
   }else if(mvv$numFit==2){             
             if(numcateg == 2){
                x = cbind(rowCoords[,1],rowCoords[,2])
                plot2CategLine(mvv,x,AtLeastR2,line=linesVoronoi,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar,PchVar=PchVar,LabValVar=LabValVar)
             }else if(QuitNotPredicted == TRUE){
                       x = cbind(rowCoords[,1],rowCoords[,2])
                       
                       plot2CategLine(mvv,x,AtLeastR2,line=linesVoronoi,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar,PchVar=PchVar,LabValVar=LabValVar)                 
                   }else{
                      plot.voronoiprob(mvv,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar,PchVar=PchVar,AtLeastR2=AtLeastR2,lines=linesVoronoi,LabValVar=LabValVar)
                   }
           }else{
              plot.voronoiprob(mvv,LabelVar=LabelVar,CexVar=CexVar,ColorVar=ColorVar,PchVar=PchVar,AtLeastR2=AtLeastR2,lines=linesVoronoi,LabValVar=LabValVar)
           }
}
