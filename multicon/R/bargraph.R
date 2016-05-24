bargraph <-
function(DV, grp=NULL, barFUN=mean, errFUN=c("ci", "se", "sd"), sides=2, conf=.95, ...) {
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

    # For a DV with no grouping variables (single bar graph)
  if(is.null(grp)) {
    DV <- na.omit(DV)
    res <- do.call(barFUN, list(DV))
    grph <- barplot(res, ...)
    if(!is.null(errFUN)) {
      e <- do.call(errFUN, list(DV))
      e <- ifelse(res <= 0, -e, e)
      if(sides==1) {
        arrows(grph, res, grph, res+e, angle=90, code=2, length=.08)
      }
      if(sides==2) {
        arrows(grph, res+e, grph, res-e, angle=90, code=3, length=.08)
      }
    }
  }
    # For a DV with one grouping variable (multiple bars - one-way ANOVA)
  if(!is.null(grp) & !is.list(grp)) {
    dat <- data.frame(DV, grp)
    dat <- dat[complete.cases(dat),]
    res <- tapply(dat[,1], dat[,2], barFUN)
    grph <- barplot(res, ...)
    if(!is.null(errFUN)) {
      e <- tapply(dat[,1], dat[,2], errFUN)
      e <- ifelse(res <= 0, -e, e)
      if(sides==1) {
        arrows(grph, res, grph, res+e, angle=90, code=2, length=.08)
      }
      if(sides==2) {
        arrows(grph, res+e, grph, res-e, angle=90, code=3, length=.08)
      }
    }
  }
    # For a DV with multiple grouping variables (multiple bars - multi-way ANOVA)
  if(is.list(grp)) {
    if(length(unique(unlist(lapply(grp, length)))) != 1) {stop("Grouping variables must be the same length.")}
    if(length(DV) != lapply(grp, length)[[1]]) {stop("DV must be the same length as the grouping variables.")}
    dat <- cbind(DV, matrix(unlist(grp), nrow=length(DV), byrow=F))
    if(sum(is.na(dat)) > 0) {stop("Please remove missing values in DV and IV first.")}
    res <- as.vector(tapply(DV, grp, barFUN))
    grph <- barplot(res, ...)
    if(!is.null(errFUN)) {
      e <- as.vector(tapply(DV, grp, errFUN))
      e <- ifelse(res <= 0, -e, e)
      if(sides==1) {
        arrows(grph, res, grph, res+e, angle=90, code=2, length=.08)
      }
      if(sides==2) {
        arrows(grph, res+e, grph, res-e, angle=90, code=3, length=.08)
      }
    }
  }
}
