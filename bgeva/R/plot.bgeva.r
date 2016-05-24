plot.bgeva <- function(x, ...){
   
    ss.plot <- x$gam.fit
    ind <- 1:length(ss.plot$coefficients)                         
    ss.plot$coefficients <- x$coefficients[ind]
    ss.plot$Vp <- x$Vb[ind,ind]
    ss.plot$edf <- diag(x$F)[ind]

    #ss.plot$smooth[[select]]$term
    
    plot.gam(ss.plot, ...)
          
}
   
