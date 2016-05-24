###Write GWR results
##Author: Binbin Lu

writeGWR<-function(x,fn="GWRresults")
{
   if(class(x) != "gwrm") stop("It's not a gwm object")
   fn1<-paste(fn,".txt",sep="")
   #fn2<-paste(fn,".csv",sep="")
   sink(fn1)
   print.gwrm(x)
   sink()
  # writeGWR.csv(x, fn=fn2)
   writeGWR.shp(x,fn=fn)
   invisible(x)
}
writeGWR.shp<-function(x,fn="GWRresults")
{
   if(class(x) != "gwrm") stop("It's not a gwm object")
   SDF<-x$SDF
   if (is(SDF, "SpatialPointsDataFrame"))
     writePointsShape(SDF,fn=fn, max_nchar= 256)
   else if (is(SDF, "SpatialPolygonsDataFrame"))
     writePolyShape(SDF, fn=fn,max_nchar= 256)
   invisible(SDF)
}