plot.Mtabs <-
function(x,...,method="normal") {
     o.Mtabs<- x
     veg<- o.Mtabs$veg
     y<- o.Mtabs$transf.r
     if(method == "normal") {
     plottab(veg,o.Mtabs$order.rel,o.Mtabs$order.sp,o.Mtabs$order.relgr,o.Mtabs$order.spgr,y)
     }
     else {
     plottabl(veg,o.Mtabs$order.rel,o.Mtabs$order.sp,o.Mtabs$order.relgr,o.Mtabs$order.spgr,y)
     }
     }
