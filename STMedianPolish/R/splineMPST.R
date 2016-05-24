#' Median polish Spline.
#'
#' The "splineMPST" is dessigned to represent the variability of effects of spatio - temporal data on a surface, from robust median polish algoritm and planar interpolation.
#' @usage splineMPST(Grid,Ef_t,MPST,eps, maxiter)
#' @param Grid grid with the coordinates in space "x", "y", "z", where will be viewed trend.
#' @param Ef_t it's the temporal scenary to look trend.
#' @param MPST object of class \code{\link{ConstructMPst}}.
#' @param eps real number greater than \code{0}, default 0.01. A tolerance for convergence.
#' @param maxiter the maximum number of iterations, default 10.
#' @return Data frame, where columns show the trend in each spatio - temporal location.
#' @references Berke, O. (2001). \emph{Modified median polish kriging and its application to the wolfcamp - aquifer data.} Environmetrics, 12(8):731-748.\href{http://onlinelibrary.wiley.com/doi/10.1002/env.495/abstract}{[link]}
#' @references Cressie, N. (1993). \emph{Statistics for spatial data.} Wiley series in probability and statistics.\href{http://www.wiley.com/WileyCDA/WileyTitle/productCd-1119115183.html}{[link]}
#' @examples 
#' ## Not run:
#' data(Metadb)
#' x<-matrix(0,1,37)
#' for(i in 1:37){
#'  x[,i] <- 2007 + (seq(0, 36)/12)[i]
#' }
#' x<-as.Date (as.yearmon(x), frac = 1)
#' time = as.POSIXct(x, tz = "GMT")
#' length(time)
#'
#' MPST<-ConstructMPst(Metadb[,-c(1:4)],time,pts=Metadb[,2:4],Delta=c(7,6,5))
#' 
#' MpSTData<-MedianPolishM(MPST,eps=0, maxiter=5, na.rm=TRUE)
#'
#' data(DemMeta)
#' xy = SpatialPoints(Metadb[,2:4],CRS(proj4string(DemMeta)))
#'
#' data(HZRMeta)
#' proj4string(HZRMeta)<-CRS(proj4string(DemMeta))
#'
#' polygon1 = polygons(HZRMeta)
#' Gridxy<- spsample(polygon1, cellsize=2000, n=300,"regular")
#'
#' Grid<-data.frame(Gridxy,over(Gridxy,DemMeta))
#' colnames(Grid)<-c("East", "North","height")
#'
#' TendenciaGrilla<-splineMPST(Grid,Ef_t=time[10:15],MPST,eps=0.01, maxiter=2)
#' 
#' IDs = paste("ID",1:length(TendenciaGrilla[,5]))
#' mydata = data.frame(values = TendenciaGrilla[,5], ID=IDs)

#' wind.ST1 = STFDF(SpatialPixels(Gridxy),time[10:15],mydata)
#' stplot(wind.ST1,col.regions=bpy.colors(40),par.strip.text = list(cex=0.7)
#'       ,main="Spline median polish: Monthly Precipitation")
#' ## End(Not run)
#' @importFrom spacetime STFDF stplot
#' @importFrom maptools readShapePoly
#' @importFrom zoo as.Date
#' @importFrom zoo as.yearmon  
#' @export 

splineMPST <-
function(Grid,Ef_t,MPST,eps=0.01, maxiter=10L){

if (dim(Grid)[2]!=3) 
stop("The interpolate grid must have 3 dimensions (x,y,z)")
if (class(MPST)!="ConstructMPst") 
stop("The MpSTData must be a a class MedianPolishM")
if (length(na.omit(Grid[,3]))!=length(Grid[,3])) 
  stop("The are NAs in Grid")

MpSTData<-MedianPolishM.default(MPST$Value, eps, maxiter, na.rm=TRUE)
pts<-MPST$pts

Delta<-MPST$Delta
ET<-rep(Ef_t,rep(dim(Grid)[1],length(Ef_t)))
Gridx<-rep(Grid[,1],length(Ef_t))
Gridy<-rep(Grid[,2],length(Ef_t))
Gridz<-rep(Grid[,3],length(Ef_t))
Gridxyz<-data.frame(Gridx,Gridy,Gridz)
colnames(Gridxyz)<-colnames(Grid)
Grid2<-data.frame(Grid,Tendencia=0)

        Grid_T<-data.frame(ET,Gridxyz,Tendencia=0)

 VDelta<-list()
        Min<-c(1:length(Delta))
        Max<-c(1:length(Delta))
        for(i in 1:length(Delta)){
        Min[i]<-min(pts[,i])}
        for(i in 1:length(Delta)){
        Max[i]<-max(pts[,i])}

#Grilla minimos y maximos del cubo

GridReal<-rbind(as.matrix(Grid[,]),as.matrix(pts[,]))
  VDeltaGrid<-list()
        MinGrid<-c(1:length(Delta))
        MaxGrid<-c(1:length(Delta))
        for(i in 1:length(Delta)){
        MinGrid[i]<-min(GridReal[,i])}
        for(i in 1:length(Delta)){
        MaxGrid[i]<-max(GridReal[,i])}

        for(i in 1:length(Delta)){
        VDelta[1:length(Delta)][[i]]<-(seq(Min[i],Max[i],abs(Min[i]-Max[i])/(Delta[i]-1)))}

        VDelta2<-VDelta

        #Remplazo de Minimos
        for(j in 1:length(Delta)){
        VDelta2[[j]][1] <- MinGrid[j]-0.1}
   
        #Remplazo de Maximos
        for(j in 1:length(Delta)){
        VDelta2[[j]][ length(VDelta2[[j]]) ] <- MaxGrid[j]+0.1}

pb <- txtProgressBar(min = 0, max = Delta[1]-1, style = 3)

for(i in 1:(Delta[1]-1)){
       for(j in 1:(Delta[2]-1)){
    for(k in 1:(Delta[3]-1)){
   sub_indexi <- which(Grid2[,1] >= VDelta2[][[1]][i] & Grid2[,1] <= VDelta2[][[1]][i+1] )
   sub_indexj <- which(Grid2[,2] >= VDelta2[][[2]][j] & Grid2[,2] <= VDelta2[][[2]][j+1] )
   sub_indexk <- which(Grid2[,3] >= VDelta2[][[3]][k] & Grid2[,3] <= VDelta2[][[3]][k+1] )
   inter1<-intersect(sub_indexj,sub_indexi)
           inter2<-intersect(inter1,sub_indexk)

   if(length(inter2)==0 ){
   }else
        {
 for(l in 1:length(inter2)){
 EX<-MpSTData$effects[[1]][i] + (Grid2[,1][inter2[l]]- VDelta[][[1]][i] )/ ( VDelta[][[1]][i+1] - VDelta[][[1]][i] )*(MpSTData$effects[[1]][i+1]-MpSTData$effects[[1]][i])
 EY<-MpSTData$effects[[2]][j] + (Grid2[,2][inter2[l]]- VDelta[][[2]][j] )/ ( VDelta[][[2]][j+1] - VDelta[][[2]][j] )*(MpSTData$effects[[2]][j+1]-MpSTData$effects[[2]][j])
 EZ<-MpSTData$effects[[3]][k] + (Grid2[,3][inter2[l]]- VDelta[][[3]][k] )/ ( VDelta[][[3]][k+1] - VDelta[][[3]][k] )*(MpSTData$effects[[3]][k+1]-MpSTData$effects[[3]][k])
                         Grid2[,4][inter2[l]]<-EX + EY + EZ  + MpSTData$overall
}
 }
    }
        }
setTxtProgressBar(pb, i)
}
close(pb)
time_c<-data.frame(MPST$time,MpSTData$effects[[4]])

for(tim in unique(ET) )
{
sub_in <- which(Grid_T[,1] == tim)
Grid_T[sub_in,5]<-Grid2[,4] + time_c[which(time_c[,1]==tim),2]
}

return(Grid_T)
}
