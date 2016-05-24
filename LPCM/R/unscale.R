
unscale<-function(x, ...){
  UseMethod("unscale")
}  



unscale.ms <- function(x,...){
 if (!x$scaled){warning("The ms-object was not fitted with scaled data, so it cannot be unscaled!")}

      cluster.center    <- sweep(x$cluster.center,2, x$scaled.by, "*")
      h                 <- x$h *x$scaled.by
      data              <- sweep(x$data,2, x$scaled.by , "*")


   return(list("cluster.center"=cluster.center, "h"=h, "data"=data)) 
}



unscale.lpc <-function(x,...){

     if (class(x)=="lpc"){ lpcobject<-x
                           splineobject<-NULL
                                }
     if (class(x)=="lpc.spline"){
                                lpcobject   <- x$lpcobject
                                splineobject<-x
                              }

     if (!lpcobject$scaled){warning("The lpcobject was not fitted with scaled data, so it cannot be unscaled!")}
      
      #if (missing(lpcobject)){lpcobject <- splineobject$lpcobject}
      
      LPC            <- sweep(lpcobject$LPC,2, lpcobject$Misc$scaled.by, "*")
      start          <- sweep(lpcobject$starting.points,2, lpcobject$Misc$scaled.by, "*")
      data           <- sweep(lpcobject$data,2, lpcobject$Misc$scaled.by , "*")
      h              <-  lpcobject$h * lpcobject$Misc$scaled.by 
      knots.coords   <- list(NULL)
      closest.coords <- list (NULL)

      if (!is.null(splineobject)){
        lk  <- length(splineobject$knots.coords)
        for (j in 1:lk){  
              knots.coords[[j]] <-  sweep(splineobject$knots.coords[[j]],1, lpcobject$Misc$scaled.by, "*")
    }
       }

      if (!is.null(splineobject) && splineobject$closest.coords!="none" ){ 
            closest.coords <- sweep(splineobject$closest.coords,2, lpcobject$Misc$scaled.by , "*")
      }
      
    return(list("LPC"=LPC, "h"=h, "data"=data, "starting.points"=start, "knots.coords"=knots.coords,"closest.coords"=closest.coords))    
 }


unscale.lpc.spline <- unscale.lpc 
