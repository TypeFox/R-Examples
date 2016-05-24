`plot.migpd` <-
function(x, main=c("Probability plot","Quantile plot","Return level plot","Histogram and density"), xlab=rep(NULL,4), nsim=1000, alpha=.05, ... ){
  
  if( class(x) != "migpd" ){
     stop("you need to use an object created by migpd")
  }
  
  if ( !missing( main ) ){
      if ( length( main ) != 1 & length( main ) != 4 ){
        stop( "main should have length 1 or 4" )
      } else if ( length( main ) == 1 ){ 
        main <- rep( main, 4 ) 
      }
  }
    
  for( i in 1:length(x$models) ) {
    plot(x$model[[i]], main= paste(rep(names(x$model[i]),4),main), xlab=xlab,nsim=nsim,alpha=0.05,...)
  }
  invisible()
}


