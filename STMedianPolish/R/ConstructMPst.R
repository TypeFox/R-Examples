
#' Construct Spatio - temporal regular data. 
#'
#' Create an  spatio - temporal object with regular data to order employing median polish technique.
#' @usage ConstructMPst(valuest,time,pts,Delta)
#' 
#' @param valuest it's a data.frame in which different columns refer to different locations, and each row reflects a particular observation time.
#' @param time indicate the time of valuest, the intervals of time must be regular.
#' @param pts it's a data.frame that hold three dimensions spatial coordinates  \code{x},  \code{y} and  \code{z}.
#' @param Delta vector with number of divisions of each spatial direction. c(Delta  \code{x}, Delta  \code{y}, Delta  \code{z}).
#' @return An object of class ConstructMPst with the following list of components:
#' @return \item{results}{average value of the set of stations into unity spatio - temporal defined for delta.}
#' @return \item{Value}{array with the results organized in dimensions defined in Delta.}
#' @return \item{valuest}{valuest}
#' @return \item{pts}{pts}
#' @return \item{time}{time}
#' @return \item{Delta}{Delta}
#' @details Table composed for coordinates of center and average of position of stations for unit tridimensional array in space - time, in which show a average value of site.
#' @references Berke, O. (2001). \emph{Modified median polish kriging and its application to the wolfcamp - aquifer data.} Environmetrics, 12(8):731-748.\href{http://onlinelibrary.wiley.com/doi/10.1002/env.495/abstract}{[link]}
#' @examples 
#' ## Not run:
#' data(Metadb)
#' x<-matrix(0,1,37)
#' for(i in 1:37){
#'  x[,i] <- 2007 + (seq(0, 36)/12)[i]
#' }
#' x<-as.Date (as.yearmon(x), frac = 1)
#' time = as.POSIXct(x, tz = "GMT")
#' 
#' MPST<-ConstructMPst(Metadb[,-c(1:4)],time,pts=Metadb[,2:4],Delta=c(7,6,5))
#' ## End(Not run)
#' @importFrom reshape2 melt
#' @importFrom zoo as.Date
#' @importFrom zoo as.yearmon 
#' @export
#' 
ConstructMPst <-
  function(valuest,time,pts,Delta)
  {
    if (!is.data.frame(valuest)) 
      stop("data is not a data.frame")
    if (length(Delta)!=3) 
      stop("Delta must have the division of the space (x,y,z)")
    if (dim(pts)[2]!=3) 
      stop("The points of the space must have coordinates (x,y,z)")
    if (dim(valuest)[1]!=dim(pts)[1])
      stop("The valuest must have the same number of rows as coordinates")
    if (length(na.omit(pts[,1]))!=length(pts[,1])||length(na.omit(pts[,2]))!=length(pts[,2]) || length(na.omit(pts[,3]))!=length(pts[,3])) 
      stop("The are NAs in pts")
        
    datast<-data.frame(time=rep(time,rep(dim(pts)[1],length(time))),pts,values=as.matrix(as.numeric(as.matrix(valuest))))
    MpData<-list()
    Data<-rep(0,1,(Delta[1]*Delta[2]*Delta[3]))
    dim(Data)<-Delta
    Mx<-Data
    My<-Mx
    Mz<-Mx
    Mx1<-Mx
    My1<-Mx
    Mz1<-Mx
    
    VDelta<-list()
    Min<-c(1:length(Delta))
    Max<-c(1:length(Delta))
    for(i in 1:length(Delta)){
      Min[i]<-min(pts[,i])  -(max(pts[,i])-min(pts[,i]))/((Delta[i]-1)*2)}
    for(i in 1:length(Delta)){
      Max[i]<-max(pts[,i])  +(max(pts[,i])-min(pts[,i]))/((Delta[i]-1)*2)}
    
    Mvalue<-matrix(0,length(unique(datast[,1])),length(Data))
    MCoorx<-Mvalue
    MCoory<-Mvalue
    MCoorz<-Mvalue
    MCoorx1<-Mvalue
    MCoory1<-Mvalue
    MCoorz1<-Mvalue
    q=0
    
    for(i in 1:length(Delta)){
      VDelta[1:length(Delta)][[i]]<-(seq(Min[i],Max[i],abs(Min[i]-Max[i])/(Delta[i])))}
    
    pb <- txtProgressBar(min = 0, max = length(unique(datast[,1])), style = 3)
    
    for(timei in unique(datast[,1]))
    {
      q=q+1
      sub_i <- which(datast[,1] == timei)
      
      for(i in 1:Delta[1]){
        for(j in 1:Delta[2]){
          for(k in 1:Delta[3]){
            sub_indexi <- which(datast[sub_i,2] >= VDelta[][[1]][i] & datast[sub_i,2] <= VDelta[][[1]][i+1])
            sub_indexj <- which(datast[sub_i,3] >= VDelta[][[2]][j] & datast[sub_i,3] <= VDelta[][[2]][j+1])
            sub_indexk <- which(datast[sub_i,4] >= VDelta[][[3]][k] & datast[sub_i,4] <= VDelta[][[3]][k+1])
            
            inter1<-intersect(sub_indexj,sub_indexi)
            inter2<-intersect(inter1,sub_indexk)
            
            if(length(inter2)==0 )
            {
              Data[i,j,k]<-NA
              Mx[i,j,k]<-mean(c(VDelta[][[1]][i],VDelta[][[1]][i+1]))
              My[i,j,k]<-mean(c(VDelta[][[2]][j],VDelta[][[2]][j+1]))
              Mz[i,j,k]<-mean(c(VDelta[][[3]][k],VDelta[][[3]][k+1]))
              Mx1[i,j,k]<-Mx[i,j,k]
              My1[i,j,k]<-My[i,j,k]
              Mz1[i,j,k]<-Mz[i,j,k]
              
            }else if(length(which(datast[sub_i,5][inter2]!='NA'))==0)
            {
              Data[i,j,k]<-NA
              Mx[i,j,k]<-sum(datast[sub_i,2][inter2],na.rm = TRUE)/length(which(datast[sub_i,2][inter2]!='NA'))
              My[i,j,k]<-sum(datast[sub_i,3][inter2],na.rm = TRUE)/length(which(datast[sub_i,3][inter2]!='NA'))
              Mz[i,j,k]<-sum(datast[sub_i,4][inter2],na.rm = TRUE)/length(which(datast[sub_i,4][inter2]!='NA'))
              Mx1[i,j,k]<-mean(c(VDelta[][[1]][i],VDelta[][[1]][i+1]))
              My1[i,j,k]<-mean(c(VDelta[][[2]][j],VDelta[][[2]][j+1]))
              Mz1[i,j,k]<-mean(c(VDelta[][[3]][k],VDelta[][[3]][k+1]))
              
            }else
            {
              
              #==========================================================
              #===========WID DEM and Values=============================
              #==========================================================
              
              MC<-datast[sub_i,5][inter2]!='NA'
              Mx[i,j,k]<-sum(datast[sub_i,2][inter2*MC],na.rm = TRUE)/length(which(datast[sub_i,2][inter2*MC]!='NA'))
              My[i,j,k]<-sum(datast[sub_i,3][inter2*MC],na.rm = TRUE)/length(which(datast[sub_i,3][inter2*MC]!='NA'))
              Mz[i,j,k]<-sum(datast[sub_i,4][inter2*MC],na.rm = TRUE)/length(which(datast[sub_i,4][inter2*MC]!='NA'))
              Mx1[i,j,k]<-mean(c(VDelta[][[1]][i],VDelta[][[1]][i+1]))
              My1[i,j,k]<-mean(c(VDelta[][[2]][j],VDelta[][[2]][j+1]))
              Mz1[i,j,k]<-mean(c(VDelta[][[3]][k],VDelta[][[3]][k+1]))
              
              if(dim((na.omit(datast[sub_i,4:5][inter2,])))[1]!=1 )
              {
                distsValue = sqrt(apply(( t(as.matrix(na.omit(datast[sub_i,2:4][inter2*MC,])))-c(Mx[i,j,k],My[i,j,k],Mz[i,j,k]))^2, 2,sum))
                distsWeigthAll<-sum((distsValue^-2.5))
                xy1<-datast[sub_i,5][inter2*MC]
                Data[i,j,k]<-sum((distsValue^-2.5)*xy1[!is.na(xy1)]/distsWeigthAll)
              }else
              {
                Data[i,j,k]<-sum(datast[sub_i,5][inter2*MC],na.rm=TRUE)
              }
            }
          }
        }
      }
      
      Mvalue[which(unique(datast[,1])==timei),]<-melt(Data)[,4]
      MCoorx[which(unique(datast[,1])==timei),]<-melt(Mx)[,4]
      MCoory[which(unique(datast[,1])==timei),]<-melt(My)[,4]
      MCoorz[which(unique(datast[,1])==timei),]<-melt(Mz)[,4]
      MCoorx1[which(unique(datast[,1])==timei),]<-melt(Mx1)[,4]
      MCoory1[which(unique(datast[,1])==timei),]<-melt(My1)[,4]
      MCoorz1[which(unique(datast[,1])==timei),]<-melt(Mz1)[,4]
      setTxtProgressBar(pb, q)
    }
    close(pb)
    
    Date<-rep(unique(datast[,1]),rep(prod(Delta),length(unique(datast[,1]))))
    
    MpData$results<-data.frame(Date,
                               XCenter=as.numeric(t(MCoorx1)), 
                               YCenter=as.numeric(t(MCoory1)),
                               ZCenter=as.numeric(t(MCoorz1)),
                               Value=as.numeric(t(Mvalue)),
                               Xmean=as.numeric(t(MCoorx)), 
                               Ymean=as.numeric(t(MCoory)),
                               Zmean=as.numeric(t(MCoorz)))
    MpData$Value<-array(MpData$results[,5],c(Delta,length(time)))
    MpData$valuest<-valuest
    MpData$pts<-pts
    MpData$time<-time
    MpData$Delta<-Delta
    
    class(MpData) <- "ConstructMPst"
    return(MpData)
}


