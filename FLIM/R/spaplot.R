spaplot <-
function(x, response, grouping=NULL, ylim=NULL, col=NULL, ...){
  dat <- x$df
  ids <- unique(dat[, 1])
  times <- x$times
  measure <- rep(0, length(times))
  if(is.null(ylim)){
    minO <- min(dat[, response])
    maxO <- max(dat[, response])
  } else {
    minO <- ylim[1]
    maxO <- ylim[2]
  }
  plot(times, measure, ylim=c(minO, maxO), type="n", axes=T, ...)
  dat.obs <- dat[dat$obs.type==1, ]
  if(is.null(grouping)){
    if(length(col)==1){col <- c(col, col)}
    if(is.null(col)){col <- c("black", "black")}
    for(i in 1:length(ids)){
      lines(times, dat[dat[, 1]==ids[i], response], lty=2, col=col[2], ...)
      lines(dat.obs[dat.obs[, 1]==ids[i], 2],
            dat.obs[dat.obs[, 1]==ids[i], response], lty =1, col=col[1], ...)
    }
    cat(paste("Spaghetti plots created for variable:", response))
    cat("\nFull drawn: observed values.")
    cat("\nStipulated: hypothetical values.\n")
    cat("\n")
  } else {
    if(is.null(col)){col <- 1:50}
    groups <- unique(dat[, grouping])
    # Ranking of nr of ids in each group level. For added visibility the biggest
    # group. Must match to group as table() rearanges alphanumerically.
    obsranking <- sort(rank(table(dat[, grouping]), ties.method="first"),
                       decreasing=TRUE)
    plotseq <- match(names(obsranking), groups)
    for(k in obsranking) {
      kdat <- dat[dat[, grouping]==groups[k], ]
      kdat.obs <- kdat[kdat$obs.type==1, ]
      kids <- unique(kdat[, 1])
      for(i in 1:length(kids)){
        lines(times, kdat[kdat[, 1]==kids[i], response], lty=2, col=col[k], ...)
        lines(kdat.obs[kdat.obs[, 1]==kids[i], 2],
              kdat.obs[kdat.obs[, 1]==kids[i], response], lty=1, col=col[k], ...)
      }
    }
    cat(paste("Spaghetti plots created for variable:", response))
    cat("\nFull drawn: observed values.")
    cat("\nStipulated: hypothetical values.\n")
    cat(paste("\nResponse variable separated by:"), grouping, "which has",
        length(groups), "levels.\n")
    cat(paste0( paste("Group:", groups), " has color ",
                col[1:length(groups)],"."))
    cat("\n")
  }
}
