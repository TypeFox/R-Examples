plot.SemiParSampleSel <- function(x, eq, ...){
 
 if(eq==1 && x$X1.d2==1) stop("There is no model component to plot.")   
 if(eq==2 && x$X2.d2==1) stop("There is no model component to plot.")   
 if(eq==3 && x$X3.d2==1) stop("There is no model component to plot.")    
 if(eq==4 && x$X4.d2==1) stop("There is no model component to plot.")   
 if(eq==5 && x$X5.d2==1) stop("There is no model component to plot.")   

 if(eq==1){ ss.plot <- x$gam1
            ind <- 1:x$X1.d2 
            } 
 if(eq==2){ ss.plot <- x$gam2
            ind <- (x$X1.d2+1):(x$X1.d2+x$X2.d2) 
            }
           
 if(eq==3){ ss.plot <- x$gam3
            ind <- (x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2) } 

 if(eq==4){ ss.plot <- x$gam4
            ind <- (x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2) }
            
 if(eq==5){ ss.plot <- x$gam5
            ind <- (x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2) }            
                                      
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$sig2 <- 1
           ss.plot$edf <- diag(x$F)[ind]
   

  plot.gam(ss.plot, ...) 
       
}



   
   
   
   
   
   
   
   