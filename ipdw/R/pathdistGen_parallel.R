#'@name pathdistGen
#'@title Generate a stack of path distance raster objects
#'@author Joseph Stachelek
#'@param spdf SpatialPointsDataFrame object
#'@param costras RasterLayer cost raster
#'@param range numeric. Range of interpolation neighborhood
#'@param step numeric. Number of sub loops to manage memory during raster processing.
#'@param yearmon character. String specifying the name of the spdf
#'@return RasterStack object of path distances
#'@import raster
#'@import gdistance
#'@export
#'@examples
#'spdf<-data.frame(rnorm(2))
#'xy<-data.frame(x=c(4,2),y=c(8,4))
#'coordinates(spdf)<-xy
#'
#'m<-matrix(NA,10,10)
#'costras<-raster(m,xmn=0,xmx=ncol(m),ymn=0,ymx=nrow(m))
#'costras[]<-runif(ncell(costras),min=1,max=10)
#'#introduce spatial gradient
#'for(i in 1:nrow(costras)){
#'costras[i,]<-costras[i,]+i
#'costras[,i]<-costras[,i]+i
#'}
#'
#'rstack<-pathdistGen(spdf,costras,100)



# 'pathdistGen2'<-function(spdf,costras,range,step=16,yearmon="default",paralleltf=TRUE){
#   
#   ipdw.range<-range/res(costras)[1]/2 #this is a per cell distance
#   #start interpolation#####
#   #calculate conductances hence 1/max(x)
#   trans<-transition(costras,function(x)1/max(x),directions=16)
#   i=1
#   coord<-spdf[i,]
#   A<-accCost(trans,coord)
#   dist<-hist(A,plot=F)$breaks[2]
#   if(dist<ipdw.range){
#     dist=ipdw.range    
#   }
# 
# if(paralleltf==TRUE){
#   library(parallel)
#   A.func<-function(trans,coord){
#     A<-accCost(trans,coord)
#     A2<-reclassify(A,c(dist,+Inf,NA,ipdw.range,dist,ipdw.range)) #raster cells are 1 map unit
#     A3<-((ipdw.range/A2)^2)
#     reclassify(A3,c(-Inf,1,0))
#   }
#   #coerce spdf to list
#   spdf.list<-list()
#   for(i in 1:length(spdf)){
#     spdf.list[[i]]<-spdf[i,]
#   }
#   cores<-detectCores()
#   cl<-makeCluster(cores-2)
#   a<-Sys.time()
#   clusterExport(cl,c("spdf.list","trans","accCost","reclassify","dist","accCost","coordinates","cellFromXY","ncell","xmin","xmax","ymin","ymax","projection","setValues","A.func","ipdw.range"))
#   A4.stack<-parLapply(cl,spdf.list,function(i)A.func(trans,i))
#   stopCluster(cl)
#   a-Sys.time()
# }
# 
#   rstack<-stack(A4.stack)
#   rstack<-reclassify(rstack,cbind(Inf,NA))
#     
#   return(rstack)
# }