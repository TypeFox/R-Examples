#
# vim:set ff=unix expandtab ts=2 sw=2:
# This test checks the code published along with Sierra et al. (2012, Geos. Model Devel. 5: 1045)
test.SoilRPaper1=function(){
     require(RUnit)
     
     attr(ICBMModel,"ex") #Shows the code
     attr(ICBMModel,"ex")() #Runs the example and shows Figure 2
     
     #Figure 3
     attr(TwopFeedbackModel, "ex") #Shows the code
     attr(TwopFeedbackModel, "ex")() #Runs the example and shows Figure 3
     
     #Figure 4. Implementing the RothC model from scratch
     attr(RothCModel, "ex") #Shows the code
     attr(RothCModel, "ex")() #Runs the example and shows Figure 4
     
     #Figure 5
     #Note: Since this part requires a large dataset, it is not included yet in the test.
#      library(sp)
#      library(raster)
#      library(ncdf)
#      library(colorspace)
#      mypal=rev(sequential_hcl(10))
#      
#      setwd(".../Data/") #Set to appropriate directory where accompanying data is stored
#      
#      PETb=(brick("PET_PT.WATCH.MonthlyMean.1980.2001.nc"))
#      Pb=(brick("Precip.WATCH.MonthlyMean.1950.2001.nc"))
#      Tbair=(brick("Tair.WATCH.MonthlyMean.1950.2001.nc"))
#      
#      PET=as.array(PETb)
#      P=as.array(Pb)
#      Tair=as.array(Tbair)
#      
#      TCent1=fT.Century1(Temp=Tair-273)
#      WCent=fW.Century(PPT=P,PET=PET)
#      
#      CDIm=TCent1*WCent
#      CDI=brick(CDIm)
#      
#      CDIa=calc(CDI,mean)
#      
#      
#      plot(CDIa,ylab="Latitude",xlab="Longitude",col=mypal,xaxt="n",yaxt="n")
     
}
