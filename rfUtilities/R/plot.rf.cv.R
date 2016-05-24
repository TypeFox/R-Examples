#' @title Plot random forests cross-validation
#' @description Plot function for rf.cv object 
#'
#' @param  x        A rf.cv object
#' @param  type     Which result to evaluate c("cv","model")
#' @param  stat     Which statistic to plot c("users.accuracy", "producers.accuracy", "kappa", "oob")  
#' @param  ...      Additional arguments passed to plot
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#'                    
#' @method plot rf.cv 
#'
#' @export    	     
plot.rf.cv <- function(x, type = "cv", stat = "producers.accuracy", ...) {
  plot.class <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- 1:nrow(x)
    dots[["y"]] <- sort(x[,1])
    dots[["type"]] <- "n"
    dots[["ylim"]] <- c(min(x),max(x)) 
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  "index"
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- stat
  	  do.call("plot", dots)
  	    graphics::axis(side=1, at = as.numeric(as.factor(rownames(x$effect.size))), 
  	                   labels = rownames(x$effect.size))
          for(i in 1:ncol(x)) { graphics::lines(dots[["x"]], sort(x[,i]), col=i) }
  	        graphics::legend("bottomright", legend=colnames(x), col=1:nrow(x), 
  	                         lty = rep(1,nrow(x)), bg="white")   
  }
    
  plot.kappa <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- stats::smooth.spline(1:nrow(x), sort(x[,"kappa"]))$x
    dots[["y"]] <- stats::smooth.spline(1:nrow(x), sort(x[,"kappa"]))$y
    dots[["type"]] <- "n"
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  ""
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "Kappa"
  	  do.call("plot", dots)
        graphics::lines(dots[["x"]], dots[["y"]]) 
  }

  plot.oob <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- 1:nrow(x)
    dots[["y"]] <- sort(x[,"OOB"])
    dots[["type"]] <- "n"
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  ""
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "OOB Error"
  	  do.call("plot", dots)
        for(i in 1:(ncol(x)-1)) { 
		  graphics::lines(1:nrow(x), sort(x[,i]), col=i)
        }
      graphics::legend("bottomright", legend = colnames(x)[1:(ncol(x)-1)],  
  	                   col=1:nrow(x), lty = rep(1,ncol(x)-1), bg="white")   		
  }
  
    if(type == "cv" & stat == "users.accuracy") 
      { dat <- x$cross.validation$cv.users.accuracy 
    } else if(type == "cv" & stat == "producers.accuracy")
      {  dat <- x$cross.validation$cv.producers.accuracy 
    } else if(type == "cv" & stat == "kappa" | stat == "oob")
      {  dat <- x$cross.validation$cv.oob 	
    } else if(type == "model" & stat == "users.accuracy")
      {  dat <- x$model$model.users.accuracy   
    } else if(type == "model" & stat == "producers.accuracy")
      {  dat <- x$model$model.producers.accuracy     
   }  else if(type == "model" & stat == "kappa" | stat == "oob")
      {  dat <- x$model$model.oob } 
	  
    if( stat == "users.accuracy" | stat == "producers.accuracy" ) {
      plot.class( dat, ...)
    } else if( stat == "oob") {
      plot.oob(dat, ...)
    } else if( stat == "kappa") {
      plot.kappa(dat, ...)
    }
} # end
