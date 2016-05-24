predict.SemiParSampleSel  <- function(object, eq, ...){
                    
 if(eq==1){ ss.pred <- object$gam1
            ind <- 1:object$X1.d2 
            } 
 if(eq==2){ ss.pred <- object$gam2
            ind <- (object$X1.d2+1):(object$X1.d2+object$X2.d2) 
            }
 if(eq==3){ ss.pred <- object$gam3
            ind <- (object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2) }       
            
 if(eq==4){ ss.pred <- object$gam4
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2) }   
            
 if(eq==5){ ss.pred <- object$gam5
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2) }               
                                                           
           ss.pred$coefficients <- object$coefficients[ind]
           ss.pred$Vp <- object$Vb[ind,ind] 
           ss.plot$sig2 <- 1

  predict.gam(ss.pred, ...)

}

