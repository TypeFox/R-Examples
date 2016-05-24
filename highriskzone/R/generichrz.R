#' Print Brief Details of a high-risk zone
#' 
#' Prints a very brief description of a high-risk zone.
#'
#' A very brief description of the highriskzone x is printed.
#' This is a method for the generic function \code{\link[base]{print}}.
#' 
#' @param x high-risk zone (object of class "\code{highriskzone}")
#' @param ... ignored
#' @method print highriskzone
#' @export 
#' @seealso \code{\link[base]{print}}, \code{\link{summary.highriskzone}}

print.highriskzone <- function(x,...){
  
  cat("high-risk zone of type", x$typehrz, " \n")
  cat("criterion:", x$criterion, " \n")
  cat("cutoff:", x$cutoff, " \n")
  
}



#' Summary of a high-risk zone
#' 
#' Prints a useful summary of a high-risk zone.
#' 
#' A useful description of the highriskzone object is printed.
#' This is a method for the generic function \code{\link[base]{summary}}.
#'
#' @param object high-risk zone (object of class "\code{highriskzone}")
#' @param ... ignored
#' @method summary highriskzone
#' @export 
#' @seealso \code{\link[base]{summary}}, \code{\link{print.highriskzone}}

summary.highriskzone <- function(object, ...){
  
  hrz <- object$zone
  hrz$m[is.na(hrz$m)] <- FALSE
  
  cat("high-risk zone of type", object$typehrz, " \n")
  cat("criterion:", object$criterion, " \n")
  cat("cutoff:", object$cutoff, " \n \n")
  
  cat("threshold:", object$threshold, " \n")
  if(!is.na(object$calccutoff)){cat("calculated cutoff:", object$calccutoff, " \n")}
  if(!is.null(object$covmatrix)){
    cat("estimated covariance matrix of Gaussian kernel:", object$covmatrix[1,], " \n")
    cat("                                               ", object$covmatrix[2,], " \n \n")
  }
  
  cat("area of the high-risk zone:", area.owin(hrz), " \n")
  
}


#' Plot a high-risk zone
#' 
#' Plot a high-risk zone.
#' 
#' This is the plot method for the class \code{highriskzone}.
#'
#' @param x high-risk zone (object of class "\code{highriskzone}")
#' @param ... extra arguments passed to the generic \code{\link[graphics]{plot}} function
#' @param pattern spatial point pattern for which the highriskzone was determined.
#' @param win observation winodw
#' @param plotpattern logical flag; if \code{TRUE}, the point pattern is plotted.
#' @param plotwindow logical flag; if \code{TRUE}, the observation window is plotted.
#' @param windowcol the color used to plot the observation window
#' @param usegpclib logical flag; if \code{TRUE}, the observation window is transformed in a 
#' polygonal window (object of class "\code{owin}" and of type "\code{polygonal}").
#' See \code{\link[spatstat]{as.polygonal}}
#' @param zonecol the colour used to plot the high-risk zone.
#' @method plot highriskzone
#' @export 
#' @seealso \code{\link[graphics]{plot}}, for examples see \code{\link{det_hrz}}

plot.highriskzone <- function(x, ..., pattern=NULL, win=NULL, plotpattern=FALSE, plotwindow=FALSE, windowcol="white", 
                              usegpclib=FALSE, zonecol="grey"){
  
  if(usegpclib){
    hrz <- x$zone
    hrz$m[is.na(hrz$m)] <- FALSE
    zone <- as.polygonal(hrz)
    plot(zone, col=c(windowcol, zonecol), ...)
  }else{
    plot(x$zone, col=c(windowcol, zonecol), ...)
  }
  if(plotwindow){
    if(!is.null(win)){
      plot(win, add=TRUE)
    }else{
      if(!is.null(pattern)){
        plot(pattern$window, add = TRUE)
      }else{stop("pattern or win must be given to plot observation window.")}
    }
  }
  if(plotpattern){
    if(!is.null(pattern)){
      plot(pattern, add = TRUE)
    }else{stop("pattern must be given to plot it.")}
  }
  
}
