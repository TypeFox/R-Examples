#'  Plot Categorical Longitudinal Data
#'
#'  Plot the empirical distribution of categorical longitudinal data.
#'
#' See http://simulx.webpopix.org/mlxr/catplotmlx/ for more details.      
#' @param r a data frame with a column \samp{id}, a column \samp{time}, 
#' a column with values and possibly a column \samp{group}.
#' @param breaks one of:
#' \itemize{
#'   \item a vector giving the breakpoints,
#'   \item a single number giving the number of segments.
#' }
#' @examples
#' \dontrun{
#'   catModel <- inlineModel("
#'   [LONGITUDINAL]
#'   input =  {a,b}
#'   EQUATION:
#'   lp1=a-b*t
#'   lp2=a-b*t/2
#'   DEFINITION:
#'   y = {type=categorical, categories={1,2,3}, 
#'   logit(P(y<=1))=lp1, logit(P(y<=2))=lp2}
#'   ")
#'   
#'   p1  <- c(a=8,b=0.2)
#'   y  <- list(name='y', time=seq(0, 100, by=4))   
#'   g1 <- list(size=100, parameter=p1)
#'   res <- simulx(model=catModel, output=y, group=g1)
#'   plot1 <- catplotmlx(res$y)
#'   print(plot1)
#'   
#'   plot2 <- catplotmlx(res$y,breaks=seq(-2,102,by=8)) 
#'   print(plot2)
#'   
#'   plot3 <- catplotmlx(res$y,breaks=5) 
#'   print(plot3)
#'   
#'   p2  <- c(a=6,b=0.3)
#'   g2 <- list(size=100, parameter=p2)
#'   res <- simulx(model=catModel, output=y, group=list(g1,g2))
#'   plot4 <- catplotmlx(res$y) 
#'   if( require("gridExtra") ){
#'     grid.arrange(plot4[[1]],plot4[[2]],ncol=2)
#'   }
#' }
#' @importFrom ggplot2 ggplot aes geom_polygon xlab ylab ylim ggtitle scale_fill_manual
#' @importFrom graphics hist
#' @importFrom grDevices hsv
#' @export         
catplotmlx <- function(r, breaks=NULL)
{
  r.name <- attr(r,"name")
  names(r)[names(r)==r.name] <- "y"
  if (is.null(breaks)){
    t <- sort(unique(r$time))
    nt <- length(t)
    zt <- c(t[1]-1,t,t[nt]+1)
    zt <- (zt[1:(nt+1)] + zt[2:(nt+2)])/2
  }else{
    if (length(breaks)>1){
    zt <- breaks
    nt <- length(zt)-1
    }else{
      nt  <- breaks
      zt1 <- min(r$time)
      zt2 <- max(r$time)
      dzt <- (zt2-zt1)/(10*nt)
      zt  <- seq(zt1-dzt,zt2+dzt,length.out=(nt+1))
    }
    t <- (zt[1:(nt)] + zt[2:(nt+1)])/2
  }
    
  y <- sort(unique(r$y))
  ny <- length(y)
  z <- c(y[1]-1,y,y[ny]+1)
  br <- (z[1:(ny+1)] + z[2:(ny+2)])/2
  v <- rep(as.factor(y), each=2*nt)
  x <- rep(c(t,rev(t)),ny)
  
  color=hsv(.6,.8,.5,seq(0.3,0.9,length.out=ny))
  sfm = scale_fill_manual(name=r.name,values=color)
  
  
  if (any( "group" %in% names(r) )){
    g=as.numeric(levels(r$group))[r$group]
  }else{
    g=rep(1,length(r$id))
  }
  ng=max(g)
  p <- vector("list", ng)
  for (kg in seq(1,ng)){
    H <- matrix(nrow=nt,ncol=ny)
    for (j in (1:nt)){
      rk<-r[g==kg,]
      yj <- rk$y[which(r$time>zt[j] & r$time<zt[j+1] )]
      hj <- hist(yj,plot=FALSE,breaks=br)
      H[j,] <- hj$density
    }
    H <- cbind(0,H)
    H <- apply(H, 1, cumsum)
    
    pr <- NULL
    for (j in (1:ny)){
      pr <- c(pr,H[j,],rev(H[j+1,]))
    }
    datapoly <- data.frame(x,pr,v)    
    pk<-ggplotmlx(datapoly, aes(x=x, y=pr)) + geom_polygon(aes(fill=v, group=v)) +
      xlab("time")+ylab("probability")+ylim(c(0,1)) 
    if (ng>1)
      pk <- pk + ggtitle(paste0("group = ",kg))
    p[[kg]] <- pk +sfm
  }
  if (ng==1)
    p <- p[[1]]
  return(p)
}  

