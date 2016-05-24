##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License 
## as published by the Free Software Foundation; either version 2 
## of the License, or (at your option) any later version.
##  
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
## Public License for more details.
##   
## You may also obtain a copy of the GNU General Public License from 
## the Free Software Foundation by visiting their Web site or by writing to
##   
##
##  Free Software Foundation, Inc.
##  59 Temple Place - Suite 330
##  Boston, MA 02111-1307
##  USA
##
################################################################################


## nTrain = Number of subjects in the training set (nGr1+nGr2)
## nRep= Number of technical replications

## sdW= sqrt(experimental or technical variation)
## sdB= sqrt(biological variation)
## rho= common Pearson correlation coefficient between biomarkers 
## foldMin= Minimum value of fold changes
## sigma= Standard deviation of the normal distribution (before truncation)
##        where fold changes are generated from  
## baseExpr = A vector of length nBiom to be used as \mu_g (log2 scale, <16) 


optimiseBiomarker<-function(error,errorTol=0.05,method="RF",nTrain=100,
                            sdB=1.5,sdW=1.0,foldAvg=2.88,nRep=3
                            ##,rho=0.0
                            )
{

  ## Set ranges for the variables in argument list
  rng.errorTol<-c(0,0.5)
  rng.nTrain<-c(10, 250)
  rng.sdB<-c(0.5, 2.5)
  rng.sdW<-c(0.1, 1.5)
  rng.foldAvg<-c(1.74,6.33)
  rng.nRep<-c(1,10)
  ##rng.rho<-c(0,0.95)

  ## Check the validity of arguments 
  if (!(method %in% c("RF","SVM","LDA","KNN")))
    stop('method should be one of "RF","SVM","LDA" and "KNN"')

  if (errorTol< rng.errorTol[1] ||errorTol> rng.errorTol[2])
  stop(paste("errorTol should be between",rng.errorTol[1],"and",rng.errorTol[2]))
  
  if (nTrain< rng.nTrain[1] ||nTrain> rng.nTrain[2])
  stop(paste("nTrain should be between",rng.nTrain[1],"and",rng.nTrain[2]))

  if (sdB< rng.sdB[1] ||sdB> rng.sdB[2])
  stop(paste("sdB should be between",rng.sdB[1],"and",rng.sdB[2]))

  if (sdW< rng.sdW[1] ||sdW> rng.sdW[2])
  stop(paste("sdW should be between",rng.sdW[1],"and",rng.sdW[2]))

  if (foldAvg< rng.foldAvg[1] ||foldAvg> rng.foldAvg[2])
  stop(paste("foldAvg should be between",rng.foldAvg[1],"and",rng.foldAvg[2]))

  if (nRep< rng.nRep[1] ||nRep> rng.nRep[2])
  stop(paste("nRep should be between",rng.nRep[1],"and",rng.nRep[2]))

  ##if (rho< rng.rho[1] ||rho> rng.rho[2])
  ##stop(paste("rho should be between",rng.rho[1],"and",rng.rho[2]))

  if(length(rgl.ids()[,1])!=0)rgl.pop(id=rgl.ids()[,1])

  panel<-rp.control(title="Determining optimal number of biomarkers",
                      error=error, method=method, nTrain=nTrain,sdB=sdB,
                      sdW=sdW,foldAvg=foldAvg,nRep=nRep,##rho=rho, 
                      errorTol=errorTol,aschar=FALSE, size=c(530,180))


  ## Radio group for methods
  rp.radiogroup(panel, method,c("RF","SVM",## "LDA",
                                "KNN"),title="Method", pos=c(10,10, 75,110),action=plot3dFun)

  ## Sliders for setting error tolerance level
  rp.slider(panel, errorTol,0,0.5, title="Error level",showvalue=TRUE,
            resolution=0.01, initval=0.05, pos=c(10,85+35,150,55),action=plot3dFun)

  ## Sliders for setting other factors
  rp.slider(panel, nTrain, 10,250, resolution=10,initval=100,title="Training set size", 
            showvalue=TRUE,pos=c(200,10,150,55),action=plot3dFun)
  rp.slider(panel, sdB, 0.5,2.5,resolution=0.1,initval=2.5, title="Biological variation",
            showvalue=TRUE,pos=c(200,10+55,150,55),action=plot3dFun)
  rp.slider(panel, sdW, 0.1,1.5,resolution=0.1 ,initval=1.5, title="Experimental variation", 
            showvalue=TRUE,pos=c(200,10+55+55,150,55),action=plot3dFun)
  rp.slider(panel, foldAvg, 1.74,6.33, resolution=0.01,initval=2.88,title="Average fold change", 
            showvalue=TRUE,pos=c(200+150+20,10,150,55),action=plot3dFun)
  rp.slider(panel, nRep, 1,10, resolution=1, initval=3,title="Replication", 
            showvalue=TRUE, pos=c(200+150+20,10+55,150,55),action=plot3dFun) 
  ##rp.slider(panel, rho, 0,0, resolution=0.05, initval=0,title="Correlation", 
  ##          showvalue=TRUE,pos=c(200+150+20,10+55+55,150,55),action=plot3dFun)
}




