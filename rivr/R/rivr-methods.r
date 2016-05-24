#' @export 
summary.rivr = function(object, ...){
  if(attr(object, "simtype") == "gvf")
    summarize_gvf(object, ...)
  else if(attr(object, "simtype") == "usf")
    summarize_usf(object, ...)
  else
    stop("Attribute 'simtype' not recognized. Supported simulations are ",
         "'usf' (unsteady flow) and 'gvf' (gradually-varied flow).")
}
         
summarize_gvf = function(x, ...){
  is = unique(round(seq(1, nrow(x), length.out = 5)))
  r = data.frame(x = x$x[is], wse = x$y[is] + x$z[is])
  r["pfn"] = round((x$y[is]/attr(x, "modspec")$normal.depth - 1)*100, 2)
  r["Fr"] = x$Fr[is]
  colnames(r) = c("Distance from control section", "Water-surface elevation",
    "Percent from normal", "Froude number")
  print(r, rownames = FALSE)
  invisible(r)
}

summarize_usf = function(x, ...){
  r = list()
  sx = x[x$monitor.type == "node",]
  r[[1]] = data.frame(distance.downstream = unique(sx$distance),
    peak.flow = tapply(sx$flow, sx$node, max),
    time.to.peak = sx$time[tapply(sx$flow, sx$node, which.max)])
  colnames(r[[1]]) = gsub("\\.", " ", names(r))
  st = x[x$monitor.type == "timestep",]
  r[[2]] = data.frame(time.since.start = unique(st$time), 
    peak.flow = tapply(st$flow, st$time, max),
    distance.travelled = st$distance[tapply(st$flow, st$time, which.max)])
  colnames(r[[2]]) = gsub("\\.", " ", names(r))
  names(r) = paste("Summary by", c("timestep", "node"))
  lapply(r, print, row.names = FALSE)
  invisible(r)
}

#' @importFrom utils head
#' @export 
head.rivr = function(x, ...){
  class(x) = "data.frame"
  NextMethod(x, ...)
}

#' @importFrom utils tail
#' @export 
tail.rivr = function(x, ...){
  class(x) = "data.frame"
  NextMethod(x, ...)
}

#' @importFrom utils str
#' @export 
print.rivr = function(x, ...){
  cat("Simulation type:\n  ")
  if(attr(x, "simtype") == "gvf")
    cat("Gradually-varied flow")
  else 
    cat("Unsteady flow", paste0("(", attr(x, "engine"), ")"))
  cat("\n\nCall:\n  ", paste(deparse(attr(x, "call")), sep = "\n", 
    collapse = "\n"), "\n\n", sep = "")
  cat("\nChannel geometry:\n")
  str(attr(x, "channel.geometry"), comp.str = " ", no.list = TRUE, 
    give.head = FALSE)
  cat("\nModel specification:\n")
  str(attr(x, "modspec"), comp.str = " ", no.list = TRUE, 
    give.head = FALSE)
  cat("\nData:\n")
  str(unclass(x), give.attr = FALSE, comp.str = " ", no.list = TRUE)  
  invisible(NULL)
}

#' @importFrom graphics plot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @export 
plot.rivr = function(x, ...){
  if (attr(x, "simtype") == "gvf") {
    with(x, plot(x, y + z, type = "l", 
      xlab = "Distance from control section", ylab = "Water-surface elevation", 
      main = "Water-surface elevation profile"))
  } else {
    par(mfrow = c(1, 2))
    with(x[x$monitor.type == "node",], {
      xvals <- tapply(time, distance, function(x) return(x))
      yvals <- tapply(flow, distance, function(x) return(x))
      plot(min(unlist(xvals)):max(unlist(xvals)), 
        ylim = (c(min(unlist(yvals)), max(unlist(yvals)))),type = "n", 
        xlab = "Time since start", ylab = "Flow", 
        main = "Flow at monitored nodes")
      mapply(lines, xvals, yvals, lty = 1:length(unique(node)))
      legend("topright", paste("x =", unique(distance)), 
        lty = 1:length(unique(node)), bty = "n")
    })
    with(x[x$monitor.type == "timestep",], {
      xvals <- tapply(distance, time, function(x) return(x))
      yvals <- tapply(flow, time, function(x) return(x))
      plot(min(unlist(xvals)):max(unlist(xvals)), 
        ylim = (c(min(unlist(yvals)), max(unlist(yvals)))),type = "n", 
        xlab = "Distance downstream", ylab = "Flow", 
        main = "Flow at monitored timesteps")
      mapply(lines, xvals, yvals, lty = 1:length(unique(step)))
      legend("topright", paste("t =", unique(time)), 
        lty = 1:length(unique(step)), bty = "n")
    })
  }  
  par(mfrow = c(1, 1))
  invisible(NULL)
}
