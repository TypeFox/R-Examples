egraph <-
function(DV, grp=NULL, plotFUN=mean, errFUN=c("ci", "se", "sd"), sides=2, conf=.95, xpoints=NULL, grp.names=NULL, tick=FALSE, ylim=NULL, len=0, ...) {
    # A function for standard error of the mean
  se <- function(x) {
    x <- na.omit(x)
    res <- sqrt(var(x)/length(x))
    res
  }
    # A function for confidence intervals for the mean
  ci <- function(x) {
    x <- na.omit(x)
    alpha <- 1-(1-conf)/2
    res <- qt(alpha, length(x)-2, lower.tail=T)*se(x)
    res
  }

  if(!is.null(errFUN)) {errFUN <- match.arg(errFUN)}
  
    # For a DV with no grouping variables (single error graph)
  if(is.null(grp)) {
    DV <- na.omit(DV)
    res <- do.call(plotFUN, list(DV))
    if(is.null(errFUN)) {plot(res, pch=19, xaxt="n", ylim=ylim, ...)}
    if(!is.null(errFUN)) {
      e <- do.call(errFUN, list(DV))
      e <- ifelse(res <= 0, -e, e)
      if(is.null(ylim)) {lims <- c(min(res+e, res-e)-abs(.4*min(res+e, res-e)), max(res+e, res-e)+.4*abs(max(res+e, res-e)))}
      else(lims <- ylim)
      plot(res, ylim=lims, pch=19, xaxt="n", ...)
      if(sides==1) {
        arrows(1, res, 1, res+e, angle=90, code=2, length=len)
      }
      if(sides==2) {
        arrows(1, res+e, 1, res-e, angle=90, code=3, length=len)
      }
    }
    if(is.null(grp.names)) {grp.names <- ""}
    axis(1, at=1, labels=grp.names, tick=tick)
  }
    # For a DV with one grouping variable (multiple lines - one-way ANOVA)
  if(!is.null(grp) & !is.list(grp)) {
    dat <- data.frame(DV, grp)
    dat <- dat[complete.cases(dat),]
    res <- tapply(dat[,1], dat[,2], plotFUN)
    if(is.null(xpoints)) {
      places <- 1:length(res)
    }
    else places <- xpoints
    
    if(is.null(errFUN)) {plot(res ~ places, pch=19, xaxt="n", xlim=c(.4,.4+places[length(places)]), ylim=ylim, ...)}
    if(!is.null(errFUN)) {
      e <- tapply(dat[,1], dat[,2], errFUN)
      e <- ifelse(res <= 0, -e, e)
      if(is.null(ylim)) {lims <- c(min(res+e, res-e)-abs(.4*min(res+e, res-e)), max(res+e, res-e)+.4*abs(max(res+e, res-e)))}
      else(lims <- ylim)
      plot(res ~ places, pch=19, xaxt="n", xlim=c(.4,.4+places[length(places)]), ylim=lims, ...)
      if(sides==1) {
        arrows(places, res, places, res+e, angle=90, code=2, length=len)
      }
      if(sides==2) {
        arrows(places, res+e, places, res-e, angle=90, code=3, length=len)
      }
    }
    if(is.null(grp.names)) {grp.names <- 1:length(places)}
    axis(1, at=places, labels=grp.names, tick=tick)
  } 
    # For a DV with multiple grouping variables (multiple lines - multi-way ANOVA)
  if(is.list(grp)) {
    if(length(unique(unlist(lapply(grp, length)))) != 1) {stop("Grouping variables must be the same length.")}
    if(length(DV) != lapply(grp, length)[[1]]) {stop("DV must be the same length as the grouping variables.")}
    dat <- data.frame(DV, matrix(unlist(grp), nrow=length(DV), byrow=F))
    if(sum(is.na(dat)) > 0) {stop("Please remove missing values in DV and IV first.")}
    res <- as.vector(tapply(DV, grp, plotFUN))
    if(is.null(xpoints)) {
      places <- 1:length(res)
    }
    else places <- xpoints
    
    if(is.null(errFUN)) {plot(res ~ places, pch=19, xaxt="n", xlim=c(.4,.4+places[length(places)]), ylim=ylim, ...)}  
    if(!is.null(errFUN)) {
      e <- as.vector(tapply(DV, grp, errFUN))
      e <- ifelse(res <= 0, -e, e)
      if(is.null(ylim)) {lims <- c(min(res+e, res-e)-abs(.4*min(res+e, res-e)), max(res+e, res-e)+.4*abs(max(res+e, res-e)))}
      else(lims <- ylim)
      plot(res ~ places, pch=19, xaxt="n", xlim=c(.4,.4+places[length(places)]), ylim=lims, ...)
      if(sides==1) {
        arrows(places, res, places, res+e, angle=90, code=2, length=len)
      }
      if(sides==2) {
        arrows(places, res+e, places, res-e, angle=90, code=3, length=len)
      }
    }
    if(is.null(grp.names)) {grp.names <- 1:length(places)}
    axis(1, at=places, labels=grp.names, tick=tick)
  }  
}
