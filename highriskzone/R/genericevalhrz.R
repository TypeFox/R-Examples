


#' Print Brief Details of an evaluation of a high-risk zone
#' 
#' Prints a very brief description of the evaluation of a high-risk zone.
#'
#' A very brief description of the evaluation x of a high-risk zone is printed.
#' This is a method for the generic function \code{\link[base]{print}}.
#' 
#' @param x evaluation of a high-risk zone (object of class "\code{hrzeval}")
#' @param ... ignored
#' @method print hrzeval
#' @export
#' @seealso \code{\link[base]{print}}, \code{\link{summary.hrzeval}}

print.hrzeval <- function(x, ...){
  cat("evaluation of a high-risk zone based on", x$numberobs, "observed events \n")
  cat("number of unobserved events:", x$numberunobserved, " \n")
  cat("number of unobserved events located outside the high-risk zone:", x$numbermiss, " \n")
}


#' Summary of a the evaluation of a high-risk zone
#' 
#' Prints a useful summary of the evaluation of a high-risk zone.
#' 
#' A useful description of the hrzeval object is printed.
#' This is a method for the generic function \code{\link[base]{summary}}.
#'
#' @param object evaluation of a high-risk zone (object of class "\code{hrzeval}")
#' @param ... ignored
#' @method summary hrzeval
#' @export 
#' @seealso \code{\link[base]{summary}}, \code{\link{print.hrzeval}}

summary.hrzeval <- function(object, ...){
  cat("evaluation of a high-risk zone based on", object$numberobs, "observed events \n")
  cat("number of unobserved events:", object$numberunobserved, " \n")
  cat("number of unobserved events located outside the high-risk zone:", object$numbermiss, " \n \n")
  
  cat("fraction of unobserved events located outside the high-risk zone:", object$missingfrac, " \n")
  cat("area of the high-risk zone:", object$arearegion, " \n \n")
  
}

#' Visualize the evaluation of a high-risk zone.
#' 
#' Plot a visualization of the evaluation of a high-risk zone. At least the observation window and the unobserved events 
#' inside and outside the high-risk zone are plotted.
#' 
#' This is the plot method for the class \code{hrzeval}.
#' 
#' @param x evaluation of a high-risk zone (object of class "\code{hrzeval}")
#' @param ... extra arguments passed to the generic \code{\link[graphics]{plot}} function.
#' @param hrz (optional) high-risk zone (object of class "\code{highriskzone}")
#' @param obspp (optional) observed point pattern
#' @param plothrz logical flag; should the high-risk zone be plotted?
#' @param plotobs logical flag; should the observed point pattern be plotted?
#' @param windowcol the color used to plot the observation window
#' @param insidecol the color used to plot the unobserved events inside the high-risk zone
#' @param outsidecol the color used to plot the unobserved events outside the high-risk zone
#' @param insidepch plotting 'character' of the unobserved events inside the high-risk zone,
#'  i.e., symbol to use. This can either be a single character or an integer code for one of 
#'  a set of graphics symbols. The full set of S symbols is available with pch=0:18, see 
#'  \code{\link[graphics]{points}}.
#' @param outsidepch plotting 'character' of the unobserved events outside the high-risk zone
#' @param zonecol the color used to plot the high-risk zone
#' @param obscol the color used to plot the observed events
#' @param obspch plotting 'character' of the observed events
#' @method plot hrzeval
#' @export 
#' @seealso \code{\link[graphics]{plot}}, \code{\link{eval_hrz}}, \code{\link{plot.highriskzone}}


plot.hrzeval <- function(x, ..., hrz=NULL, obspp=NULL, plothrz=FALSE, plotobs=FALSE, windowcol="white", 
                         insidecol="blue", outsidecol="red", insidepch = 20, outsidepch = 19, 
                         zonecol="grey", obscol="black", obspch = 1){
  
  if(plotobs & is.null(obspp)){stop("observed pattern (obspp) must be given to plot it.")}
  if(plothrz & is.null(hrz)){stop("high-risk zone (hrz) must be given to plot it.")}
  
  if(plothrz){
    plot(hrz, zonecol=zonecol, win=x$out$window, plotwindow=TRUE, windowcol=windowcol, ...)
  } else {
    plot(x$out$window, col=windowcol, ...)
  }
  if(plotobs){
    plot(obspp, add=TRUE, pch=obspch, col=obscol)
  }
  plot(x$insd, col=insidecol, pch=insidepch, add = TRUE)
  plot(x$out, col=outsidecol, pch=outsidepch, add = TRUE)
  
}
