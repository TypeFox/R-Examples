plot.SemiParBIVProbit <- function(x, eq, ...){

 if(missing(eq)) stop("You must provide the equation from which smooth terms should be plotted.")


 if(eq > x$l.flist) stop("The fitted model has a smaller number of equations.") 


 if(eq==1 && x$l.sp1==0) stop("There is no model component to plot.")   
 if(eq==2 && x$l.sp2==0) stop("There is no model component to plot.")   
 if(eq==3 && x$l.sp3==0) stop("There is no model component to plot.")    
 if(eq==4 && x$l.sp4==0) stop("There is no model component to plot.")   
 if(eq==5 && x$l.sp5==0) stop("There is no model component to plot.")   
 if(eq==6 && x$l.sp6==0) stop("There is no model component to plot.")   
 if(eq==7 && x$l.sp7==0) stop("There is no model component to plot.")   
 

 if(eq==1){ ss.plot <- x$gam1
            ind <- 1:x$X1.d2 
            } 
 if(eq==2){ ss.plot <- x$gam2
            ind <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2) 
            }
 if(eq==3){ ss.plot <- x$gam3
            ind <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) }
            
 if(eq==4){ ss.plot <- x$gam4
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) }
            
 if(eq==5){ ss.plot <- x$gam5
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2) }
            
 if(eq==6){ ss.plot <- x$gam6
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2) }   
            
 if(eq==7){ ss.plot <- x$gam7
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2) }             
                                    
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$sig2 <- 1
           ss.plot$edf <- diag(x$F)[ind]
           ss.plot$scale.estimated <- FALSE 

  plot.gam(ss.plot, ...)  
       
}

