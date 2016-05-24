
plot.pscore <- function(x,
                        par.dens=NULL,
                        par.1=NULL,
                        par.0=NULL,
                        with.legend=FALSE,
                        legend.cex=0.9,
                        legend.label=NULL,
                        main=NULL,
                        ylim=NULL,
                        xlim=NULL,
                        ...)
{
  object <- x
  
  if ( any(class(object) == "pscore") ){
    
    require( "colorspace", character.only=TRUE )
    
    if (is.null(ylim))
      ylim <- c(0, ceiling(max(density(object$pscore[object$treat==1])$y,
                               density(object$pscore[object$treat==0])$y)))
    if (is.null(xlim))
      xlim <- c(floor(min(density(object$pscore[object$treat==1])$x,
                          density(object$pscore[object$treat==0])$x)),
                ceiling(max(density(object$pscore[object$treat==1])$x,
                            density(object$pscore[object$treat==0])$x)))
    
    plot(density(object$pscore),
         col="white",
         main=main,
         ylim=ylim,
         xlim=xlim,
         ...)
    
    x1 <- do.call("density",
                  c(list(x=object$pscore[object$treat==1]),
                    par.dens))  
    x0 <- do.call("density",
                  c(list(x=object$pscore[object$treat==0]),
                    par.dens))
    
    if (is.null(par.1) & is.null(par.0))
      
      par.0 <- list(lty=2)
    
    
    do.call("lines", c(list(x=x1), par.1))
    do.call("lines", c(list(x=x0), par.0))
    
    if (with.legend){
      
      if (is.null(legend.label)){
        legend.label <- c("treated", "untreated")}
      
      legend(x=xlim[2]/4,
             y=ylim[2],
             cex=legend.cex,
             lty=c(1:2),
             legend.label,
             ncol=2)
    }
  }else{
    stop("Function 'pscore' must be at least used before.")
  }
}
  
