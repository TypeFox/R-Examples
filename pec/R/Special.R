#' Drawing bootstrapped cross-validation curves and the .632 or .632plus error
#' of models. The prediction error for an optional benchmark model can be added
#' together with bootstrapped cross-validation error and apparent errors.
#' 
#' This function is invoked and controlled by \code{plot.pec}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{Special.arg} in the call to \code{plot.pec}.
#' 
#' @param x an object of class 'pec' as returned by the \code{pec} function.
#' @param y Prediction error values.
#' @param addprederr Additional prediction errors. The options are bootstrap
#' cross-validation errors or apparent errors.
#' @param models One model also specified in \code{pec} for which the
#' \code{predErr} in \code{plot.pec} is to be drawn.
#' @param bench A benchmark model (also specified in \code{pec}) for which the
#' \code{predErr} in \code{plot.pec} is to be drawn.
#' @param benchcol Color of the benchmark curve.
#' @param times Time points at which the curves must be plotted.
#' @param maxboot Maximum number of bootstrap curves to be added. Default is
#' all.
#' @param bootcol Color of the bootstrapped curves. Default is 'gray77'.
#' @param col Color of the different error curves for \code{models}.
#' @param lty Line type of the different error curves for \code{models}.
#' @param lwd Line width of the different error curves for \code{models}.
#' @return Invisible object.
#' @seealso \code{\link{plot.pec}}
Special <-  function(x,
                     y,
                     addprederr,
                     models,
                     bench,
                     benchcol,
                     times,
                     maxboot,
                     bootcol,
                     col,
                     lty,
                     lwd){
   if(length(models) != 1)
     stop("Need to choose one and only one 'models'")
   at <- prodlim::sindex(x$time,times)

   # add the bootstrap curves 
   boot <- x$BootstrapCrossValErrMat[[models]]
  if(is.null(maxboot)) maxboot <- NCOL(boot) 
   
  nix <- lapply(1:maxboot,function(b){
    lines(times,boot[b,at,drop=TRUE],col=bootcol,type="s",lwd=2.5)
  })

   # add the predErr of models
   addw <- lines(times, y[[models]][at], col=col[1], lty=lty[1], lwd=lwd[1])
   
   # add the xtra chosen whats of models
   if(!is.null(addprederr)){
     nx <- length(addprederr)
     addx <- lapply(1:nx, function(a){
       newy <- x[[addprederr[a]]][[models]]
       lines(times, newy[at], col=col[a+1], lwd=lwd[a+1], lty=lty[a+1])
     })
   }
   
   # add the benchmark models
   if(bench!= FALSE){
     addb <- lines(times, y=y[[bench]][at], col=benchcol, lwd=lwd[1])
   }
 }

