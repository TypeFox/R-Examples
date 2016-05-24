plot.flim <-
function(x, response=NULL, grouping=NULL, ylim=NULL, col=NULL, naive=T,
                      lty=1:2, ptype="mean", ylab="Response", xlab="Times", ...)
  {
  if(is.null(response)) response <- x$info$responses[1]
  if(ptype=="mean"){
    dataset <- x$df
    times <- x$times
    if(is.null(grouping)) {
      if(length(col)==1){col <- c(col, col)}
      if(is.null(col)){col <- c("black", "black")}
      mean.responses <- flimMean(x, response, grouping)
      if(is.null(ylim)){
        minO <- min(mean.responses)
        maxO <- max(mean.responses)
      } else {
        minO <- ylim[1]
        maxO <- ylim[2]
      }
      plot(times, mean.responses[, 1], ylim=c(minO, maxO), type="n", axes=T,
           ylab=ylab, xlab=xlab, ...)
      lines(times, mean.responses[, 1], lty=lty[2], col=col[2], ...)
      points(times, mean.responses[, 1], lty=lty[2], col=col[2], ...)
      cat(paste("Mean response plot created for variable:", response))
      cat("\nFull drawn: observed mean response.")
      cat("\nStipulated: flim hypothetical mean response.\n")
      cat("\n")
      if(naive==TRUE) {
        
        lines(times, mean.responses[, 2], col=col[1], lty=lty[1], )
        points(times, mean.responses[, 2], col=col[1], lty=lty[1], ...)
      }
    }
    else {
      if(is.null(col)){col <- 1:40}
      if(is.null(ylim)){
        minO <- min(flimMean(x, response, grouping))
        maxO <- max(flimMean(x, response, grouping))
      } else {
        minO <- ylim[1]
        maxO <- ylim[2]
      }
      measure <- rep(0,length(times))
      plot(times, measure, ylim=c(minO, maxO), type="n", axes=T, ylab=ylab,
           xlab=xlab,...)
      groups <- unique(dataset[, grouping]) 
      for(l in 1:length(groups)) {
        group <- groups[l]
        dat.all <- dataset[dataset[, grouping] == group, ]
        dat.obs <- dat.all[dat.all$obs.type == 1, ]
        mr.all <- tapply(dat.all[, response], dat.all[, 2], mean)
        mr.obs <- tapply(dat.obs[, response], dat.obs[, 2], mean)
        lines(times, mr.all, lty=2, col=col[l])
        points(times, mr.all, lty=2, col=col[l], ...)
        if(naive==TRUE) {
          lines(times, mr.obs, lty=1, col=col[l])
          points(times, mr.obs, lty=1, col=col[l], ...)
        }
      }
      cat(paste("Mean response plot created for variable:", response))
      cat("\nFull drawn: observed mean response.")
      cat("\nStipulated: flim hypothetical mean response.\n")
      cat(paste("\nResponse variable separated by:"), grouping, "which has",
          length(groups), "levels.\n")
      cat(paste0( paste("Group:", groups), " has color ",
                  col[1:length(groups)],"."))
      cat("\n")
    }
  }
  if(ptype=="spa"){
    spaplot(x, response, grouping, ylim, col, ylab=ylab, xlab=xlab, ...)
  }
}
