


##############################################################
# S3 plot method for objects of class rasch.mml2
plot.rasch.mml <- function( x , items=NULL , xlim=NULL , main=NULL , 
		...){
    object <- x
    theta <- object$theta.k
	if ( is.matrix(theta) ){
	   stop("Plot function is only applicable for unidimensional models")
						}
    probs <- t(object$pjk)
    if ( is.null(items) ){ items <- 1:(nrow(probs) ) }
    I <- length(items)
    if (is.null(xlim)){ xlim <- c( min(theta) , max(theta) ) }    
    xlabplot <- expression( paste( theta ))
    ylabplot <- expression(paste(P,"(",X[i]==1, "|", theta , ")" ))    
	if (object$irtmodel == "ramsay.qm"){
			xlabplot <- expression( paste( "log " ,theta ))	
						}	
    graphics::plot( theta , as.vector(probs[1,]) , type="l" , lty=1 , xlab=xlabplot , 
        ylab=ylabplot , xlim=xlim , ylim=c(0,1) , main=main , ... )
    for (ii in 2:I){
        graphics::lines( theta , as.vector(probs[ii,]) , type="l" , lty=ii)
                    }
            }
##############################################################