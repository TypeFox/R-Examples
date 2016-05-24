#' Plot Multiple Conditional/Unconditional Indices 
#' @param con a concatenated series of data frames produced by SEAsic functions for (up to 6) conditional indices 
#' @param u a concatenated series of scalars representing (up to 6) unconditional indices, which can be entered directly or called in from previous SEAsic estimation of such indices
#' @param s scalar representing the standard deviation of equated scores in the overall population, or on any (sub)population of interest (e.g., synthetic population) (default is 1, which should be used when plotting unstandardized indices)
#' @param d unstandardized difference that matters
#' @param connames a row vector of quoted names for the conditional indices
#' @param unames a row vecotr of quoted names for the unconditional indices
#' @param ylab a quoted label for the y axis (default is "Equating Dependence")
#' @param xlab a quoted label for the x axis (default is "Score Scale")
#' @param colc a color designation for the conditional indices (default is a series of rainbow colors the length of u)
#' @param colu a color designation for the unconditional index lines (default is a series of rainbow colors matching colc)
#' @param pchc a row vector of shape designation for the conditional index data points (default is filled circle)
#' @param ltyu a row vector of line types for the unconditional index lines (default is straight line)
#' @param yleg a scalar indicating the value of the y axis at which the legend should sit (default is the maximum conditional index value for the first data frame in con)
#' @param xleg a scalar indicating the value of the x axis at which the legend should sit (default is the minimum value of the score scale in the first data frame in con)
#' @param ylim a row vector length of 2 with the minimum and maximum values for the y axis (default set at finite values in the plot)
#' @author Anne Corinne Huggins-Manley 
#' @seealso \code{\link{pc}}
#' @examples
#' #Obtaining and plotting the unstandardized RSD(x) and RESD for all five subpopulations in the example data set, ex.data
#' rsd_g1 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],d=.5)
#' rsd_g2 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],d=.5)
#' rsd_g3 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],d=.5)
#' rsd_g4 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],d=.5)
#' rsd_g5 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5)
#' resd_g1 <- resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],f=ex.data[,8])
#' resd_g2 <- resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],f=ex.data[,8])
#' resd_g3 <- resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],f=ex.data[,8])
#' resd_g4 <- resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],f=ex.data[,8])
#' resd_g5 <- resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],f=ex.data[,8])
#' 
#' pcu(con=c(rsd_g1,rsd_g2,rsd_g3,rsd_g4,rsd_g5),u=c(resd_g1,resd_g2,resd_g3,resd_g4,resd_g5),d=.5,connames=c("RSD Group 1","RSD Group 2","RSD Group 3","RSD Group 4","RSD Group 5"),unames=c("RESD Group 1","RESD Group 2","RESD Group 3","RESD Group 4","RESD Group 5"),ylim=c(0,8),yleg=8)
#' 
#' #Obtaining and plotting the standardized RSD(x) and RESD for all five subpopulations in the example data set, ex.data
#' srsd_g1 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],d=.5,s=4.2)
#' srsd_g2 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],d=.5,s=4.2)
#' srsd_g3 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],d=.5,s=4.2)
#' srsd_g4 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],d=.5,s=4.2)
#' srsd_g5 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5,s=4.2)
#' sresd_g1<-resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],f=ex.data[,8],s=4.2)
#' sresd_g2<-resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],f=ex.data[,8],s=4.2)
#' sresd_g3<-resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],f=ex.data[,8],s=4.2)
#' sresd_g4<-resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],f=ex.data[,8],s=4.2)
#' sresd_g5<-resd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],f=ex.data[,8],s=4.2) 
#' 
#' pcu(con=c(srsd_g1,srsd_g2,srsd_g3,srsd_g4,srsd_g5),u=c(sresd_g1,sresd_g2,sresd_g3,sresd_g4,sresd_g5),d=.5,s=4.2,connames=c("RSD Group 1","RSD Group 2","RSD Group 3","RSD Group 4","RSD Group 5"),unames=c("RESD Group 1","RESD Group 2","RESD Group 3","RESD Group 4","RESD Group 5"),ylim=c(0,2),yleg=2,ylab="Standardized Equating Dependence")
#' 
#' #Obtaining and plotting the unstandardized RMSD(x) and REMSD for all five subpopulations in the example data set, ex.data
#' rmsd <- rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),w=c(.1,.2,.4,.2,.1),d=.5)
#' remsd <- remsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),f=ex.data[,8],w=c(.1,.2,.4,.2,.1))
#' 
#' pcu(con=c(rmsd),u=c(remsd),d=.5,connames=c("RMSD"),unames=c("REMSD"),ylim=c(0,4),yleg=4)
#' @note All indices must be measured on the same scale (e.g., all unstandardized indices), and all conditional indices must have an unconditional counterpart in the plot.
#' @return plot of multiple conditional and unconditional indices
#' @export

pcu<- function(con,u,s,d,connames,unames,ylab,xlab,colc,colu,pchc,ltyu,yleg,xleg,ylim){
  
  con <- as.data.frame(con)  
  
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(ylab))
    ylab <- "Equating Dependence"
  if(missing(s))
    s <- 1
  if(missing(colc))
    colc <- rainbow(n=length(u))
  if(missing(colu))
    colu <- colc
  if(missing(pchc))
    pchc <- c(rep(19,length(u)))
  if(missing(ltyu))
    ltyu <- c(rep(1,length(u)))
  if(missing(yleg))
    yleg <- max(con[,2])
  if(missing(xleg))
    xleg <- min(con[,1])
  if(missing(ylim))
    ylim <- NULL
  
  
  if(length(u)==1){
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])
                  legend(x=xleg,y=yleg,c(connames[1],unames[1],"DTM"),col=c(colc[1],colu[1],"black"),pch=c(pchc[1],-1,-1),lty=c(0,ltyu[1],5))}
  else{   
  if(length(u)==2){
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])
                  points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
                  abline(u[2],0,col=colu[2],lty=ltyu[2])
                  legend(x=xleg,y=yleg,c(connames[1],connames[2],unames[1],unames[2],"DTM"),
                         col=c(colc[1],colc[2],colu[1],colu[2],"black"),pch=c(pchc[1],pchc[2],-1,-1,-1),lty=c(0,0,ltyu[1],ltyu[2],5))}
  else{
  if(length(u)==3){
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])  
                  points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
                  abline(u[2],0,col=colu[2],lty=ltyu[2])
                  points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
                  abline(u[3],0,col=colu[3],lty=ltyu[3])
                  legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],unames[1],unames[2],unames[3],"DTM"),
                      col=c(colc[1],colc[2],colc[3],colu[1],colu[2],colu[3],"black"),pch=c(pchc[1],pchc[2],pchc[3],-1,-1,-1,-1),lty=c(0,0,0,ltyu[1],ltyu[2],ltyu[3],5))}
  else{
  if(length(u)==4){      
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])  
                  points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
                  abline(u[2],0,col=colu[2],lty=ltyu[2])
                  points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
                  abline(u[3],0,col=colu[3],lty=ltyu[3])
                  points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
                  abline(u[4],0,col=colu[4],lty=ltyu[4])  
                  legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],unames[1],unames[2],unames[3],unames[4],"DTM"),
                      col=c(colc[1],colc[2],colc[3],colc[4],colu[1],colu[2],colu[3],colu[4],"black"),
                      pch=c(pchc[1],pchc[2],pchc[3],pchc[4],-1,-1,-1,-1,-1),lty=c(0,0,0,0,ltyu[1],ltyu[2],ltyu[3],ltyu[4],5))}
  else{
  if(length(u)==5){  
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])  
                  points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
                  abline(u[2],0,col=colu[2],lty=ltyu[2])
                  points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
                  abline(u[3],0,col=colu[3],lty=ltyu[3])
                  points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
                  abline(u[4],0,col=colu[4],lty=ltyu[4])
                  points(type="b",con[,1],con[,10],col=colc[5],pch=pchc[5])
                  abline(u[5],0,col=colu[5],lty=ltyu[5])
                  legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],connames[5],unames[1],unames[2],unames[3],unames[4],unames[5],"DTM"),
                      col=c(colc[1],colc[2],colc[3],colc[4],colc[5],colu[1],colu[2],colu[3],colu[4],colu[5],"black"),
                      pch=c(pchc[1],pchc[2],pchc[3],pchc[4],pchc[5],-1,-1,-1,-1,-1,-1),lty=c(0,0,0,0,0,ltyu[1],ltyu[2],ltyu[3],ltyu[4],ltyu[5],5))}
  else{
  if(length(u)==6){  
      pmcuplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
                  abline(d/s,0,col="black",lty=2)
                  abline(u[1],0,col=colu[1],lty=ltyu[1])  
                  points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
                  abline(u[2],0,col=colu[2],lty=ltyu[2])
                  points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
                  abline(u[3],0,col=colu[3],lty=ltyu[3])
                  points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
                  abline(u[4],0,col=colu[4],lty=ltyu[4])
                  points(type="b",con[,1],con[,10],col=colc[5],pch=pchc[5])
                  abline(u[5],0,col=colu[5],lty=ltyu[5])
                  points(type="b",con[,1],con[,12],col=colc[6],pch=pchc[6])
                  abline(u[6],0,col=colu[6],lty=ltyu[6])
                  legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],connames[5],connames[6],unames[1],unames[2],unames[3],unames[4],unames[5],unames[6],"DTM"),
                      col=c(colc[1],colc[2],colc[3],colc[4],colc[5],colc[6],colu[1],colu[2],colu[3],colu[4],colu[5],colu[6],"black"),
                      pch=c(pchc[1],pchc[2],pchc[3],pchc[4],pchc[5],pchc[6],-1,-1,-1,-1,-1,-1,-1),lty=c(0,0,0,0,0,0,ltyu[1],ltyu[2],ltyu[3],ltyu[4],ltyu[5],ltyu[6],5))}
  }}}}}
  
  
  
  
  return(pmcuplot)
}