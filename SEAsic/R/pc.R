#' Plot a Single or Multiple Conditional Indices 
#' @param con a concatenated series of data frames produced by SEAsic functions for (up to 6) conditional indices 
#' @param n scalar indicating the number of conditional indices in the plot
#' @param s scalar representing the standard deviation of equated scores in the overall population, or on any (sub)population of interest (e.g., synthetic population) (default is 1, which should be used when plotting unstandardized indices)
#' @param d unstandardized difference that matters
#' @param connames a row vector of quoted names for the conditional indices
#' @param ylab a quoted label for the y axis (default is "Equating Dependence")
#' @param xlab a quoted label for the x axis (default is "Score Scale")
#' @param colc a color designation for the conditional indices (default is a series of n rainbow colors)
#' @param pchc a row vector of shape designation for the conditional index data points (default is filled circle)
#' @param yleg a scalar indicating the value of the y axis at which the legend should sit (default is the maximum conditional index value for the first data frame in con)
#' @param xleg a scalar indicating the value of the x axis at which the legend should sit (default is the minimum value of the score scale in the first data frame in con)
#' @param ylim a row vector length of 2 with the minimum and maximum values for the y axis (default set at finite values in the plot)
#' @author Anne Corinne Huggins-Manley 
#' @seealso \code{\link{pcu}}
#' @examples
#' #Obtaining and plotting the unstandardized RSD(x) for all five subpopulations in the example data set, ex.data
#' rsd_g1 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],d=.5)
#' rsd_g2 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],d=.5)
#' rsd_g3 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],d=.5)
#' rsd_g4 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],d=.5)
#' rsd_g5 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5)
#' 
#' pc(con=c(rsd_g1,rsd_g2,rsd_g3,rsd_g4,rsd_g5),n=5,d=.5,connames=c("RSD Group 1","RSD Group 2","RSD Group 3","RSD Group 4","RSD Group 5"),ylim=c(0,4),yleg=4)
#' 
#' #Obtaining and plotting the standardized RSD(x) for all five subpopulations in the example data set, ex.data
#' srsd_g1 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],d=.5,s=4.2)
#' srsd_g2 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,4],d=.5,s=4.2)
#' srsd_g3 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,5],d=.5,s=4.2)
#' srsd_g4 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,6],d=.5,s=4.2)
#' srsd_g5 <- rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5,s=4.2)
#' 
#' pc(con=c(srsd_g1,srsd_g2,srsd_g3,srsd_g4,srsd_g5),n=5,d=.5,connames=c("RSD Group 1","RSD Group 2","RSD Group 3","RSD Group 4","RSD Group 5"),s=4.2,ylim=c(0,2),yleg=2,ylab="Standardized Equating Dependence")
#' 
#' #Obtaining and plotting the unstandardized RMSD(x) for all five subpopulations in the example data set, ex.data
#' rmsd <- rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),w=c(.1,.2,.4,.2,.1),d=.5)
#' 
#' pc(con=c(rmsd),n=1,d=.5,connames=c("RMSD"),ylim=c(0,4),yleg=4,ylab="Root Mean Square Differences")

#' @note 
#' \itemize{
#' \item{All indices must be measured on the same scale (e.g., all unstandardized indices).}
#' \item{Plotting a single conditional index can also be done via functions to obtain the values of any conditional index, but there are more options to customize the plots via the pc function.}
#' }
#' @return plot of multiple conditional indices
#' @export

pc<- function(con,n,s,d,connames,ylab,xlab,colc,pchc,yleg,xleg,ylim){ 
  
  con <- as.data.frame(con)  
  
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(ylab))
    ylab <- "Equating Dependence"
  if(missing(s))
    s <- 1
  if(missing(colc))
    colc <- rainbow(n=n)
  if(missing(pchc))
    pchc <- c(rep(19,n))
  if(missing(yleg))
    yleg <- max(con[,2])
  if(missing(xleg))
    xleg <- min(con[,1])
  if(missing(ylim))
    ylim <- NULL
  
  
  
  if(n==1){
    pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
    abline(d/s,0,col="black",lty=2)
    legend(x=xleg,y=yleg,c(connames[1],"DTM"),col=c(colc[1],"black"),pch=c(pchc[1],-1),lty=c(0,5))}
  else{   
    if(n==2){
      pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
      abline(d/s,0,col="black",lty=2)
      points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
      legend(x=xleg,y=yleg,c(connames[1],connames[2],"DTM"),
             col=c(colc[1],colc[2],"black"),pch=c(pchc[1],pchc[2],-1),lty=c(0,0,5))}
    else{
      if(n==3){
        pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
        abline(d/s,0,col="black",lty=2)
        points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
        points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
        legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],"DTM"),
               col=c(colc[1],colc[2],colc[3],"black"),pch=c(pchc[1],pchc[2],pchc[3],-1),lty=c(0,0,0,5))}
      else{
        if(n==4){      
          pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
          abline(d/s,0,col="black",lty=2)
          points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
          points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
          points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
          legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],"DTM"),
                 col=c(colc[1],colc[2],colc[3],colc[4],"black"),
                 pch=c(pchc[1],pchc[2],pchc[3],pchc[4],-1),lty=c(0,0,0,0,5))}
        else{
          if(n==5){  
            pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
            abline(d/s,0,col="black",lty=2)
            points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
            points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
            points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
            points(type="b",con[,1],con[,10],col=colc[5],pch=pchc[5])
            legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],connames[5],"DTM"),
                   col=c(colc[1],colc[2],colc[3],colc[4],colc[5],"black"),
                   pch=c(pchc[1],pchc[2],pchc[3],pchc[4],pchc[5],-1),lty=c(0,0,0,0,0,5))}
          else{
            if(n==6){  
              pmcplot <- plot(type="b",x=con[,1],y=con[,2],xlab=xlab,ylab=ylab,col=colc[1],pch=pchc[1],ylim=ylim)
              abline(d/s,0,col="black",lty=2)
              points(type="b",con[,1],con[,4],col=colc[2],pch=pchc[2])
              points(type="b",con[,1],con[,6],col=colc[3],pch=pchc[3])
              points(type="b",con[,1],con[,8],col=colc[4],pch=pchc[4])
              points(type="b",con[,1],con[,10],col=colc[5],pch=pchc[5])
              points(type="b",con[,1],con[,12],col=colc[6],pch=pchc[6])
              legend(x=xleg,y=yleg,c(connames[1],connames[2],connames[3],connames[4],connames[5],connames[6],"DTM"),
                     col=c(colc[1],colc[2],colc[3],colc[4],colc[5],colc[6],"black"),
                     pch=c(pchc[1],pchc[2],pchc[3],pchc[4],pchc[5],pchc[6],-1),lty=c(0,0,0,0,0,0,5))}
          }}}}}
  
  
  
  
  return(pmcplot)
}