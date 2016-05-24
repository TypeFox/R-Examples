plot.lpc <- function(x, type, unscale=TRUE, lwd=1, datcol="grey60",   datpch=21, masscol=NULL, masspch=15, curvecol=1, splinecol=3, projectcol=4, startcol=NULL, startpch=NULL,...){

     
   object <- x
   if (class(object)=="lpc.spline"){
       splineobject <- object
       lpcobject    <- object$lpcobject
   } else {
       lpcobject <- object
       splineobject <-NULL
   }
      
   if (missing(type)){
     if (class(object)=="lpc"){type="curve"}
     else if (class(object)=="lpc.spline" && object$project==FALSE){type="spline"}
     else if (class(object)=="lpc.spline" && object$project){type=c("spline","project")}
   }                                                   
          
      
   if ("project" %in% type  && class(object)=="lpc" || "project" %in% type &&  class(object)=="lpc.spline" && object$closest.branch=="none"  ){
       splineobject <-    lpc.spline(lpcobject, project=TRUE)
   } else if ("spline" %in% type  && class(object)=="lpc"){
       splineobject <-    lpc.spline(lpcobject)
   }
        
   n        <- dim(lpcobject$data)[1]
   d        <- dim(lpcobject$data)[2]
   branch   <- as.factor(lpcobject$P[,"branch"])
   nbranch  <- nlevels(branch)
   lc       <- length(curvecol)   
   sc       <- length(splinecol)
   pc       <- length(projectcol)  
   stnames  <- dimnames(lpcobject$start)[[1]]
   stdepth  <- as.numeric(stnames)   
   maxdepth <- max(stdepth)
  
   
   if (is.null(masscol)){
     masscol<-as.numeric(c(2,4,1)[1:maxdepth])
   }
   if (is.null(startcol)){
     startcol<-as.numeric(c(6,5,"grey20")[1:maxdepth])
   }
   
   mc       <- length(masscol)
   stc      <- length(startcol)  
     
  if (!is.null(startpch)){
          if (length(stnames)!=length(startpch)){
             stnames<-rep(startpch[1], length(stnames))
           } else {
             stnames<-startpch
           }
  }
   
  if (lpcobject$scaled && unscale){
        u       <- if (is.null(splineobject)) {unscale(lpcobject)} else  {unscale(splineobject)}
        Xi      <- u$data
        fit     <- u$LPC
        start   <- u$starting.points
        if ("spline" %in% type) knots.coords <- u$knots.coords
        if ("project" %in% type) closest.coords <- u$closest.coords
  } else {
        Xi    <-  lpcobject$data
        fit   <-  lpcobject$LPC
        start <-  lpcobject$starting.points
        if ("spline" %in% type) knots.coords <- splineobject$knots.coords
        if ("project" %in% type) closest.coords <- splineobject$closest.coords 
}  
   
   
if (d==2){
   # eqscplot(Xi, col=datcol, pch=datpch,...)
   plot(Xi, col=datcol, pch=datpch,...)
   if ("curve" %in% type){
          if (maxdepth==1){     
               for (j in 0:(nbranch-1)){ 
                   lines(fit[branch==j,], lwd=lwd, col= if (lc>=nbranch) curvecol[j+1] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol)
              }
          } else {
               for (j in 0:(nbranch-1)){ 
                     lines(fit[branch==j,], lwd=lwd, col= if (lc>=maxdepth) curvecol[stdepth[j+1]] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol)
               }         
          }
   }       
   
          
  if ("mass" %in% type){
      if (mc== dim(fit)[1]){
             points(fit, col=masscol,  pch=masspch )
      } else {
         for (j in 0:(nbranch-1)){ 
            points(fit[branch==j,], pch=masspch, col=  if (mc==maxdepth)  masscol[as.numeric(dimnames(start)[[1]])[j+1] ]  else  if  (mc>=nbranch) masscol[j+1] else if
          (mc>1) rep(masscol,ceiling(nbranch/mc))[j+1] else masscol)
         }
      }   
   }
   if ("start" %in% type){
      points(start, col = if (stc == maxdepth) startcol[stdepth] else startcol, pch=stnames)
   } 
   
   if ("spline" %in% type){
       for (j in 0:(nbranch-1)){ 
           lines(knots.coords[[j+1]][1,],knots.coords[[j+1]][2,], lwd=lwd, col= if (sc>=nbranch) splinecol[j+1] else if (sc>1) rep(splinecol,ceiling(nbranch/sc))[j+1] else splinecol)
         }
   }    
   if ("project" %in% type){
       for (i in 1:n){
         x1 <-closest.coords[i,1]
         y1<- closest.coords[i,2]
         x2 <- Xi[i,1]
         y2 <- Xi[i,2]
        segments(x1,y1,x2,y2,  col= if (pc==n) projectcol[i] else if (pc>=nbranch) projectcol[splineobject$closest.branch[i]+1] else if (pc>1)  rep(projectcol,ceiling(nbranch/pc))[splineobject$closest.branch[i]+1] else projectcol  ) # 09/10/09
       }
   }
  
      
   
} else if (d==3){
  
   #require(scatterplot3d)
   plotlpc3 <- scatterplot3d::scatterplot3d(Xi, color=datcol, pch=datpch,...)
   if ("curve" %in% type){ 
         for (j in 0:(nbranch-1)){ 
             plotlpc3$points3d(fit[branch==j,], lwd=lwd, col=if (lc>=nbranch) curvecol[j+1] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol, type="l")
         }         
       }


   
   if ("mass" %in% type){
          if (mc== dim(fit)[1]){
              plotlpc3$points3d(fit, col=masscol,  pch=masspch )
          } else {  
              for (j in 0:(nbranch-1)){ 
              plotlpc3$points3d(fit[branch==j,], pch=masspch, col= if (mc==maxdepth)  masscol[as.numeric(dimnames(start)[[1]])[j+1] ]  else  if (mc>=nbranch) masscol[j+1] else if (mc>1) rep(masscol,ceiling(nbranch/mc))[j+1] else masscol, type="p")
              }
          }    
        }
   
    if ("start" %in% type){
      plotlpc3$points3d(start, col = if (stc == maxdepth) startcol[stdepth] else startcol, pch=stnames)
   }  


   
  if ("spline" %in% type){
       for (j in 0:(nbranch-1)){
           plotlpc3$points3d(t(knots.coords[[j+1]]), lwd=lwd, col= if (sc>=nbranch) splinecol[j+1] else if (sc>1) rep(splinecol,ceiling(nbranch/sc))[j+1] else splinecol, type="l")
       }
   }
  if ("project" %in% type){
  
       for (i in 1:n){
         x1 <- closest.coords[i,1]
         y1<-  closest.coords[i,2]
         z1 <- closest.coords[i,3] 
         x2 <- Xi[i,1]
         y2 <- Xi[i,2]
         z2 <- Xi[i,3]
         plotlpc3$points3d(c(x1,x2),c(y1,y2),c(z1,z2), type="l",
              col= if (pc==n) projectcol[i]  else if (pc>=nbranch) projectcol[splineobject$closest.branch[i]+1] else if (pc>1)  rep(projectcol,ceiling(nbranch/pc))[splineobject$closest.branch[i]+1] else projectcol)  # 09/10/2009
       }
   }    
 
   
   
}  else if (d<=16){
        pairs(fit,
              panel= if ("curve" %in%  type) "lines" else "points",
              labels= dimnames(lpcobject$data)[[2]],
              lwd=lwd,
              col=curvecol,
              ...
              )
  
  #X11()
  #require(lattice)
  #splom(lpcobject$LPC)
  ##splom(lpcobject$LPC)#, type=if ("curve" %in%  type) "l" else "p",labels= dimnames(lpcobject$data)[[2]] )
  
}  else {
  cat("Data set has too many dimensions to be plotted in a pairs plot.\n\n")
}


invisible()

}
  


plot.lpc.spline <- plot.lpc
