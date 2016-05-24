


plot.invariance.alignment <- function( x , ... ){
	mod3 <- x
    graphics::par( mfrow=c(2,1))
graphics::plot( 1:(mod3$Niter[1]) , mod3$fopt.history[ 1:(mod3$Niter[1]) ,1 ] , type="l" ,
        xlim = c(0,max(mod3$Niter)) , ylim= range( mod3$fopt.history[,1]  , na.rm=TRUE) ,
        xlab="Iteration" , ylab= "Optimization Function" ,
        main="Optimization History LAMBDA")
graphics::plot( 1:(mod3$Niter[2]) , mod3$fopt.history[ 1:(mod3$Niter[2]) ,2] , type="l" ,
        xlim = c(0,max(mod3$Niter)) , ylim= range( mod3$fopt.history[,2]  , na.rm=TRUE) ,
        xlab="Iteration" , ylab= "Optimization Function" ,
        main="Optimization History NU")        
    graphics::par(mfrow=c(1,1))
        }