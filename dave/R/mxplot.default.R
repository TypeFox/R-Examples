mxplot.default <-
function(veg,rmember,use,y=1,...)  {
     o.mxplot<- matrixplot(veg,rmember,use,y=1,...)
     o.mxplot$call<- match.call()
     cat("Call:\n") 
     class(o.mxplot) <- "mxplot"
     print(o.mxplot$call)
     o.mxplot
     }
